#' Static Count Distribution Model
#'
#' Static (IID) count distribution model for intermittent demand, following
#' Kolassa (2016). The method fits several candidate distributions --- Poisson,
#' hurdle-shifted Poisson, negative binomial, and hurdle-shifted negative
#' binomial --- to the observed series and selects the best by AIC. A mixture
#' option that blends all four predictive distributions is also available.
#'
#' @param formula Model specification.
#' @param distr Distribution choice: one of `"auto"`, `"pois"`, `"hsp"`,
#'   `"nbinom"`, `"hsnb"`, or `"mixture"`.
#' @param hot_start Logical. If `TRUE`, leading zeros are removed from the
#'   time series before fitting.
#' @param ... Not used.
#'
#' @references
#' Kolassa, S. (2016). Evaluating predictive count data distributions in retail
#' sales forecasting. *International Journal of Forecasting*, 32(3), 788--803.
#'
#' @return A model specification.
#'
#' @importFrom fabletools new_model_class new_specials new_model_definition
#' @importFrom tsibble measured_vars
#' @importFrom rlang abort arg_match is_integerish
#' @importFrom distributional dist_sample
#' @importFrom nloptr nloptr
#' @importFrom stats dpois dnbinom rpois rnbinom runif var setNames
#' @export
STATICDISTR <- function(formula, distr = c("auto", "pois", "hsp", "nbinom", "hsnb", "mixture"), hot_start = FALSE, ...) {
  distr <- arg_match(distr)

  staticdistr_model <- new_model_class(
    "STATICDISTR",
    train = train_staticdistr,
    specials = new_specials(
      xreg = staticdistr_no_xreg
    )
  )
  new_model_definition(staticdistr_model, {{ formula }}, distr = distr, hot_start = hot_start, ...)
}

train_staticdistr <- function(.data, specials, distr, hot_start = FALSE, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by STATICDISTR.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by STATICDISTR.")
  }

  # Identify the distributions to be fitted
  if (distr %in% c("auto", "mixture")) {
    c("nbinom", "pois", "hsnb", "hsp")
  } else {
    distr
  }

  # Apply Croston's decomposition
  decomp <- crostons_decomp(y)
  occurrence <- decomp$occurrence
  shifted_demand <- decomp$demand - 1

  fit_results <- list()
  if ("pois" %in% distributions) {
    fit_results[["pois"]] <- staticdistr_fit_pois(y)
  }
  if ("hsp" %in% distributions) {
    fit_results[["hsp"]] <- staticdistr_fit_hsp(y, occurrence, shifted_demand)
  }
  if ("nbinom" %in% distributions) {
    fit_results[["nbinom"]] <- staticdistr_fit_nbinom_result(y)
  }
  if ("hsnb" %in% distributions) {
    fit_results[["hsnb"]] <- staticdistr_fit_hsnb(y, occurrence, shifted_demand)
  }
  fit_results

  # Extrapolate the distribution via AIC or BIC
  aic <- vapply(fit_results, function(x) x$aic, numeric(1))
  bic <- vapply(fit_results, function(x) x$bic, numeric(1))
  mles <- staticdistr_collect_mles(fit_results)
  pred_distr <- if (distr == "mixture") to_eval else names(which.min(aic))
  fitted_mean <- staticdistr_expected_value(pred_distr, fit_results)
  fitted <- rep(fitted_mean, length(y))

  structure(
    list(
      aic = aic,
      bic = bic,
      mles = mles,
      fit_results = fit_results,
      pred_distr = pred_distr,
      fitted = fitted,
      residuals = y - fitted
    ),
    class = "STATICDISTR"
  )
}

#' @export
forecast.STATICDISTR <- function(object, new_data, specials = NULL, times = 10000, ...) {
  h <- nrow(new_data)
  if (!is_integerish(times) || times <= 0) {
    abort("`times` must be a positive integer.")
  }

  sim <- staticdistr_simulate(object, h = h, times = as.integer(times))
  dist_sample(as.list(as.data.frame(sim)))
}

#' @export
generate.STATICDISTR <- function(x, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  sim <- staticdistr_simulate(x, h = h, times = 1L)
  new_data$.sim <- as.numeric(sim[1, ])
  new_data
}

#' @export
fitted.STATICDISTR <- function(object, ...) {
  object$fitted
}

#' @export
residuals.STATICDISTR <- function(object, ...) {
  object$residuals
}

#' @export
model_sum.STATICDISTR <- function(x) {
  "STATICDISTR"
}

staticdistr_simulate <- function(object, h, times) {
  iid_n <- h * times
  preds <- object$pred_distr
  fit_results <- object$fit_results

  per <- floor(iid_n / length(preds))
  rem <- iid_n %% length(preds)
  counts <- rep(per, length(preds))
  if (rem > 0) {
    counts[seq_len(rem)] <- counts[seq_len(rem)] + 1L
  }

  samples <- numeric(0)
  for (i in seq_along(preds)) {
    di <- preds[[i]]
    ni <- counts[[i]]
    if (ni <= 0) {
      next
    }

    samples <- c(samples, staticdistr_simulate_from_distribution(di, ni, fit_results))
  }

  if (length(preds) > 1) {
    samples <- samples[sample.int(length(samples))]
  }

  matrix(samples, nrow = times, ncol = h, byrow = TRUE)
}

staticdistr_expected_value <- function(pred_distr, fit_results) {
  means <- numeric(0)

  for (di in pred_distr) {
    means <- c(means, staticdistr_distribution_mean(di, fit_results[[di]]$params))
  }

  mean(means)
}

staticdistr_fit_pois <- function(y) {
  params <- c(lambda = mean(y))
  loglik <- sum(dpois(y, params[["lambda"]], log = TRUE))
  staticdistr_new_fit_result("pois", params, loglik, length(y))
}

staticdistr_fit_hsp <- function(y, occurrence, shifted_demand) {
  params <- c(
    pzero = mean(1 - occurrence),
    lambda = if (length(shifted_demand) > 0) mean(shifted_demand) else 0
  )
  loglik <- sum(staticdistr_dhsp(y, params[["pzero"]], params[["lambda"]], log = TRUE))
  staticdistr_new_fit_result("hsp", params, loglik, length(y))
}

staticdistr_fit_nbinom_result <- function(y) {
  params <- staticdistr_fit_nbinom(y)
  loglik <- sum(dnbinom(y, params[["size"]], params[["prob"]], log = TRUE))
  staticdistr_new_fit_result("nbinom", params, loglik, length(y))
}

staticdistr_fit_hsnb <- function(y, occurrence, shifted_demand) {
  if (length(shifted_demand) > 0) {
    nbinom_params <- staticdistr_fit_nbinom(shifted_demand)
  } else {
    nbinom_params <- c(size = 100, prob = 1 - staticdistr_epsilon)
  }

  params <- c(
    pzero = mean(1 - occurrence),
    size = nbinom_params[["size"]],
    prob = nbinom_params[["prob"]]
  )
  loglik <- sum(staticdistr_dhsnb(y, params[["pzero"]], params[["size"]], params[["prob"]], log = TRUE))
  staticdistr_new_fit_result("hsnb", params, loglik, length(y))
}

staticdistr_new_fit_result <- function(distribution, params, loglik, n_obs) {
  n_params <- length(params)

  list(
    distribution = distribution,
    params = params,
    loglik = loglik,
    aic = -2 * loglik + 2 * n_params,
    bic = -2 * loglik + log(n_obs) * n_params,
    n_params = n_params
  )
}

staticdistr_collect_mles <- function(fit_results) {
  setNames(
    unlist(lapply(names(fit_results), function(dist_name) {
      fit <- fit_results[[dist_name]]
      setNames(unname(fit$params), paste0(dist_name, "_", names(fit$params)))
    })),
    unlist(lapply(names(fit_results), function(dist_name) {
      fit <- fit_results[[dist_name]]
      paste0(dist_name, "_", names(fit$params))
    }))
  )
}

staticdistr_distribution_mean <- function(distribution, params) {
  if (distribution == "pois") {
    return(params[["lambda"]])
  }

  if (distribution == "hsp") {
    return((1 - params[["pzero"]]) * (1 + params[["lambda"]]))
  }

  if (distribution == "nbinom") {
    return(params[["size"]] * (1 - params[["prob"]]) / params[["prob"]])
  }

  if (distribution == "hsnb") {
    return((1 - params[["pzero"]]) *
      (1 + params[["size"]] * (1 - params[["prob"]]) / params[["prob"]]))
  }

  abort(paste0("Unsupported distribution: ", distribution))
}




staticdistr_fit_nbinom <- function(y) {
  if (length(y) == 0 || all(y == 0)) {
    return(c(size = 100, prob = 1 - staticdistr_epsilon))
  }

  fit <- tryCatch(
    nloptr(
      x0 = c(max(mean(y), staticdistr_epsilon), 0.5),
      eval_f = function(x) -mean(dnbinom(y, x[1], x[2], log = TRUE)),
      lb = c(staticdistr_epsilon, staticdistr_epsilon),
      ub = c(Inf, 1 - staticdistr_epsilon),
      opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || is.null(fit$solution)) {
    mu <- mean(y)
    sigmasq <- var(y)
    if (!is.na(sigmasq) && sigmasq > mu + staticdistr_epsilon) {
      size <- (mu^2) / (sigmasq - mu)
    } else {
      size <- 100
    }
    prob <- min(size / (size + mu), 1 - staticdistr_epsilon)
    return(c(size = size, prob = prob))
  }

  c(size = fit$solution[1], prob = fit$solution[2])
}



dhsp <- function(x, pzero = 0.5, lambda = 1, log = FALSE) {
  if (log) {
    ifelse(x == 0, log(pzero), log(1 - pzero) + dpois(x - 1, lambda, log = TRUE))
  } else {
    ifelse(x == 0, pzero, (1 - pzero) * dpois(x - 1, lambda))
  }
}

phsp <- function(q, pzero = 0.5, lambda = 1) {
  ifelse(q < 0, 0, pzero + (1 - pzero) * ppois(q - 1, lambda))
}

qhsp <- function(p, pzero = 0.5, lambda = 1) {
  ifelse(p <= pzero, 0, qpois((p - pzero) / (1 - pzero), lambda) + 1)
}

rhsp <- function(n, pzero = 0.5, lambda = 0.1) {
  ifelse(runif(n) <= pzero, 0, 1 + rpois(n, lambda))
}



dhsnb <- function(x, pzero = 0.5, size = 1, prob = 0.5, log = FALSE) {
  if (log) {
    ifelse(x == 0, log(pzero), log(1 - pzero) + dnbinom(x - 1, size, prob, log = TRUE))
  } else {
    ifelse(x == 0, pzero, (1 - pzero) * dnbinom(x - 1, size, prob))
  }
}

phsnb <- function(q, pzero = 0.5, size = 1, prob = 0.5) {
  ifelse(q < 0, 0, pzero + (1 - pzero) * pnbinom(q - 1, size, prob))
}

qhsnb <- function(p, pzero = 0.5, size = 1, prob = 0.5) {
  ifelse(p <= pzero, 0, qnbinom((p - pzero) / (1 - pzero), size, prob) + 1)
}

rhsnb <- function(n, pzero = 0.5, size = 1, prob = 0.5) {
  ifelse(runif(n) <= pzero, 0, 1 + rnbinom(n, size, prob))
}


staticdistr_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by STATICDISTR.")
}

staticdistr_epsilon <- 1e-4
