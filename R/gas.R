#' @param formula Model specification.
#' @param distr Distribution choice: one of `"pois"`, `"nbinom"`, or `"bern"`.
#' @param level_type Character. One of `"RW"` or `"AR"`.
#' @param ... Not used.
#'
#' @references
#' Sarlo, R., Fernandes, C., & Borenstein, D. (2023). Lumpy and intermittent
#' retail demand forecasts with score-driven models. *European Journal of
#' Operational Research*, 307(3), 1146--1160.
#'
#' @return A model specification.
#'
#' @importFrom fabletools new_model_class new_specials new_model_definition model_sum
#' @importFrom tsibble measured_vars
#' @importFrom rlang abort is_integerish arg_match
#' @importFrom distributional dist_sample dist_poisson
#' @importFrom nloptr nloptr
#' @importFrom stats dpois rpois
GAS <- function(formula, ...) {
  level_type <- arg_match(level_type)

  gas_model <- new_model_class(
    "GAS",
    train = train_gas,
    specials = new_specials(
      xreg = gas_no_xreg
    )
  )
  new_model_definition(gas_model, {{ formula }}, distr = distr, level_type = level_type, ...)
}

train_gas <- function(.data, specials, distr, level_type, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by GAS.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by GAS.")
  }
  if (distr == "bern" && any(!y %in% c(0, 1))) {
    abort("`distr = \"bern\"` requires a binary response.")
  }

  period <- get_freq(.data)
  opt <- gas_optimize(y, period, distr, level_type)
  x <- opt$solution
  params <- gas_unpack_params(x, period, distr, level_type)

  rec <- gas_filter(
    y,
    params$psi0,
    params$phi,
    params$rho,
    params$xi0,
    params$k,
    period,
    distr,
    params$alpha
  )

  fitted <- rec$f
  residuals <- y - fitted
  last_score <- gas_score(y[length(y)], fitted[length(fitted)], distr, params$alpha)

  structure(
    list(
      distr = distr,
      level_type = level_type,
      period = period,
      psi0 = params$psi0,
      phi = params$phi,
      rho = params$rho,
      xi0 = if (period > 1) params$xi0 else NA_real_,
      k = if (period > 1) params$k else NA_real_,
      alpha = params$alpha,
      n = length(y),
      last_y = y[length(y)],
      last_f = fitted[length(fitted)],
      last_score = last_score,
      last_psi = rec$last_psi,
      last_xi = rec$last_xi,
      fitted = fitted,
      residuals = residuals
    ),
    class = "GAS"
  )
}

#' Forecast a GAS model
#'
#' Produces forecast distributions from a fitted GAS model using simulation.
#'
#' @inheritParams forecast.EMPDISTR
#' @param times The number of sample paths to use in estimating the forecast
#'   distribution.
#'
#' @export
forecast.GAS <- function(object, new_data, specials = NULL, times = 10000, ...) {
  h <- nrow(new_data)
  if (!is_integerish(times) || times <= 0) {
    abort("`times` must be a positive integer.")
  }

  f_first <- gas_first_step(object)
  dist_first <- gas_distribution(f_first, object$distr, object$alpha)

  if (h == 1) {
    return(dist_first)
  }

  sim <- gas_simulate(object, h, times)
  samples_rest <- as.list(as.data.frame(sim[, -1, drop = FALSE]))
  dist_rest <- dist_sample(samples_rest)

  c(dist_first, dist_rest)
}

#' Generate sample paths from a GAS model
#'
#' @param x A fitted `GAS` model object.
#' @inheritParams forecast.GAS
#' @export
generate.GAS <- function(x, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  sim <- gas_simulate(x, h, 1L)
  new_data$.sim <- as.numeric(sim[1, ])
  new_data
}

#' Extract fitted values from a GAS model
#'
#' @inheritParams forecast.GAS
#' @export
fitted.GAS <- function(object, ...) {
  object$fitted
}

#' Extract residuals from a GAS model
#'
#' @inheritParams forecast.GAS
#' @export
residuals.GAS <- function(object, ...) {
  object$residuals
}

#' @export
model_sum.GAS <- function(x) {
  paste0("GAS[", x$distr, ", ", x$level_type, "]")
}

gas_first_step <- function(object) {
  score <- object$last_score
  psi <- object$phi * object$last_psi + object$rho * score

  if (object$period > 1) {
    season <- (object$n %% object$period) + 1
    gamma <- object$last_xi[season] + object$k * score
    gas_linkinv(psi + gamma, object$distr)
  } else {
    gas_linkinv(psi, object$distr)
  }
}

gas_simulate <- function(object, h, times) {
  forecast_samples <- matrix(NA_real_, nrow = times, ncol = h)

  score_state <- rep(object$last_score, times)
  psi_state <- rep(object$last_psi, times)
  if (object$period > 1) {
    xi_state <- matrix(rep(object$last_xi, each = times), nrow = times)
  }

  for (i in seq_len(h)) {
    psi_state <- object$phi * psi_state + object$rho * score_state

    if (object$period > 1) {
      season <- ((object$n + i - 1) %% object$period) + 1
      xi_state[, season] <- xi_state[, season] + object$k * score_state
      xi_state[, -season] <- xi_state[, -season] - object$k / (object$period - 1) * score_state
      f_state <- gas_linkinv(psi_state + xi_state[, season], object$distr)
    } else {
      f_state <- gas_linkinv(psi_state, object$distr)
    }

    forecast_samples[, i] <- gas_sample(object$distr, f_state, object$alpha, times)
    score_state <- gas_score(forecast_samples[, i], f_state, object$distr, object$alpha)
  }

  forecast_samples
}

gas_filter <- function(y, psi0, phi, rho, xi0, k, period, distr, alpha) {
  if (distr == "nbinom") {
    gasFilterNbinom(y, psi0, phi, rho, xi0, k, period, alpha)
  } else if (distr == "bern") {
    gasFilterBern(y, psi0, phi, rho, xi0, k, period)
  } else {
    gasFilterPois(y, psi0, phi, rho, xi0, k, period)
  }
}

gas_score <- function(y, f, distr, alpha = NULL) {
  if (distr == "pois") {
    y - f
  } else if (distr == "nbinom") {
    (y - f) / (1 + f / alpha)
  } else {
    (1 - y) - f
  }
}

gas_linkinv <- function(eta, distr) {
  if (distr == "bern") {
    stats::plogis(eta)
  } else {
    exp(eta)
  }
}

gas_sample <- function(distr, f, alpha, times) {
  if (distr == "pois") {
    rpois(times, lambda = pmax(f, gas_epsilon))
  } else if (distr == "nbinom") {
    size <- if (is.null(alpha)) 1 else alpha
    rnbinom(times, mu = pmax(f, gas_epsilon), size = size)
  } else {
    rbinom(times, size = 1, prob = pmin(pmax(f, gas_epsilon), 1 - gas_epsilon))
  }
}

gas_distribution <- function(f, distr, alpha) {
  if (distr == "pois") {
    dist_poisson(f)
  } else if (distr == "nbinom") {
    size <- if (is.null(alpha)) 1 else alpha
    prob <- size / (size + f)
    dist_negative_binomial(size = size, prob = prob)
  } else {
    distributional::dist_binomial(size = 1, prob = pmin(pmax(f, gas_epsilon), 1 - gas_epsilon))
  }
}

gas_unpack_params <- function(x, period, distr, level_type) {
  idx <- 1
  psi0 <- x[idx]
  idx <- idx + 1

  if (level_type == "AR") {
    phi <- x[idx]
    idx <- idx + 1
  } else {
    phi <- 1
  }

  rho <- x[idx]
  idx <- idx + 1

  if (period > 1) {
    xi0 <- x[idx]
    idx <- idx + 1
    k <- x[idx]
    idx <- idx + 1
  } else {
    xi0 <- 0
    k <- 0
  }

  alpha <- NULL
  if (distr == "nbinom") {
    alpha <- exp(x[idx])
  }

  list(psi0 = psi0, phi = phi, rho = rho, xi0 = xi0, k = k, alpha = alpha)
}

gas_optimize <- function(y, period, distr, level_type) {
  mean_y <- mean(y)
  max_y <- max(y)
  psi_start <- if (distr == "bern") {
    stats::qlogis(min(max(mean_y, gas_epsilon), 1 - gas_epsilon))
  } else {
    log(max(mean_y, gas_epsilon))
  }

  if (level_type == "AR") {
    if (period > 1) {
      if (distr == "nbinom") {
        x0 <- c(psi_start, 0.8, 0.1, psi_start / 2, 0.1, log(max(mean_y, 1)))
        lb <- c(gas_psi_lower(distr), -1 + gas_epsilon, gas_epsilon, gas_psi_lower(distr), gas_epsilon, log(gas_epsilon))
        ub <- c(gas_psi_upper(distr, max_y), 1 - gas_epsilon, 1 - gas_epsilon, gas_psi_upper(distr, max_y), 1 - gas_epsilon, log(max(max_y, 1) + 100))
      } else {
        x0 <- c(psi_start, 0.8, 0.1, psi_start / 2, 0.1)
        lb <- c(gas_psi_lower(distr), -1 + gas_epsilon, gas_epsilon, gas_psi_lower(distr), gas_epsilon)
        ub <- c(gas_psi_upper(distr, max_y), 1 - gas_epsilon, 1 - gas_epsilon, gas_psi_upper(distr, max_y), 1 - gas_epsilon)
      }
    } else {
      if (distr == "nbinom") {
        x0 <- c(psi_start, 0.8, 0.1, log(max(mean_y, 1)))
        lb <- c(gas_psi_lower(distr), -1 + gas_epsilon, gas_epsilon, log(gas_epsilon))
        ub <- c(gas_psi_upper(distr, max_y), 1 - gas_epsilon, 1 - gas_epsilon, log(max(max_y, 1) + 100))
      } else {
        x0 <- c(psi_start, 0.8, 0.1)
        lb <- c(gas_psi_lower(distr), -1 + gas_epsilon, gas_epsilon)
        ub <- c(gas_psi_upper(distr, max_y), 1 - gas_epsilon, 1 - gas_epsilon)
      }
    }
  } else {
    if (period > 1) {
      if (distr == "nbinom") {
        x0 <- c(psi_start, 0.1, psi_start / 2, 0.1, log(max(mean_y, 1)))
        lb <- c(gas_psi_lower(distr), gas_epsilon, gas_psi_lower(distr), gas_epsilon, log(gas_epsilon))
        ub <- c(gas_psi_upper(distr, max_y), 1 - gas_epsilon, gas_psi_upper(distr, max_y), 1 - gas_epsilon, log(max(max_y, 1) + 100))
      } else {
        x0 <- c(psi_start, 0.1, psi_start / 2, 0.1)
        lb <- c(gas_psi_lower(distr), gas_epsilon, gas_psi_lower(distr), gas_epsilon)
        ub <- c(gas_psi_upper(distr, max_y), 1 - gas_epsilon, gas_psi_upper(distr, max_y), 1 - gas_epsilon)
      }
    } else {
      if (distr == "nbinom") {
        x0 <- c(psi_start, 0.1, log(max(mean_y, 1)))
        lb <- c(gas_psi_lower(distr), gas_epsilon, log(gas_epsilon))
        ub <- c(gas_psi_upper(distr, max_y), 1 - gas_epsilon, log(max(max_y, 1) + 100))
      } else {
        x0 <- c(psi_start, 0.1)
        lb <- c(gas_psi_lower(distr), gas_epsilon)
        ub <- c(gas_psi_upper(distr, max_y), 1 - gas_epsilon)
      }
    }
  }

  nloptr(
    x0 = x0,
    eval_f = function(x) gas_nll(x, y, period, distr, level_type),
    lb = lb,
    ub = ub,
    opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
  )
}

gas_nll <- function(x, y, period, distr, level_type) {
  params <- gas_unpack_params(x, period, distr, level_type)
  rec <- gas_filter(
    y,
    params$psi0,
    params$phi,
    params$rho,
    params$xi0,
    params$k,
    period,
    distr,
    params$alpha
  )

  if (distr == "pois") {
    -mean(dpois(y, lambda = pmax(rec$f, gas_epsilon), log = TRUE))
  } else if (distr == "nbinom") {
    -mean(dnbinom(y, mu = pmax(rec$f, gas_epsilon), size = params$alpha, log = TRUE))
  } else {
    -mean(dbinom(y, size = 1, prob = pmin(pmax(rec$f, gas_epsilon), 1 - gas_epsilon), log = TRUE))
  }
}

gas_psi_lower <- function(distr) {
  if (distr == "bern") {
    stats::qlogis(gas_epsilon)
  } else {
    log(gas_epsilon)
  }
}

gas_psi_upper <- function(distr, max_y) {
  if (distr == "bern") {
    stats::qlogis(1 - gas_epsilon)
  } else {
    log(max(max_y, 1) + 100)
  }
}

gas_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by GAS.")
}

gas_epsilon <- 1e-4
