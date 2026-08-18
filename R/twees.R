#' Tweedie Exponential Smoothing
#'
#' Exponential smoothing state space model for intermittent demand with a
#' Tweedie observation distribution. The conditional mean of the Tweedie
#' is governed by a (optionally damped) exponential smoothing process, and
#' the dispersion parameter is derived by a second (optionally damped) 
#' exponential smoothing process on the occurrence binary time series. 
#' The Tweedie family naturally models both zeros and large spikes 
#' via its compound Poisson-Gamma nature. The power parameter is optimised
#' to maximise the likelihood. The first-step forecast follows a Tweedie 
#' distribution, and multi-step forecasts are obtained by simulating 
#' from the model forward in time.
#'
#' @param formula Model specification.
#' @param damped Logical. If `TRUE` (default), the exponential smoothing
#'   component uses a damping parameter.
#' @param scaling Logical. If `TRUE` (default), the time series is divided by
#'   its maximum value before fitting and predictions are back-transformed.
#'   This improves numerical stability.
#' @param object A fitted model object.
#' @param ... Not used.
#'
#' @references
#'
#' Damato, S., Azzimonti, D., & Corani, G. (2025). Forecasting intermittent
#' time series with Gaussian Processes and Tweedie likelihood.
#' *International Journal of Forecasting*  (in press).
#' \doi{10.1016/j.ijforecast.2025.10.001}.
#'
#' @return A model specification.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#'
#' fc_ts <- ts |>
#'   model(TWEES(value)) |>
#'   forecast(h = "7 days")
#'
#' fc_ts |> print()
#'
#'
#' if (requireNamespace("ggtime", quietly = TRUE)) {
#'   library(ggtime)
#'   fc_ts |> autoplot(ts)
#' }
#' @importFrom fabletools new_model_class new_specials new_model_definition
#' @importFrom tsibble measured_vars
#' @importFrom rlang abort is_integerish
#' @importFrom distributional dist_sample
#' @importFrom nloptr nloptr
#' @importFrom stats median
#' @importFrom tweedieDistr dist_tweedie dtweedie rtweedie
#' @export
TWEES <- function(formula, damped = TRUE, scaling = TRUE, ...) {
  twees_model <- new_model_class(
    "TWEES",
    train = train_twees,
    specials = new_specials(
      xreg = twees_no_xreg
    )
  )
  new_model_definition(twees_model, {{ formula }}, damped = damped, scaling = scaling, ...)
}

#' @importFrom stats median
train_twees <- function(.data, specials, damped, scaling, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by TWEES.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by TWEES.")
  }
  if (!is.logical(damped)) {
    abort("`damped` must be a boolean.")
  }

  # Optionally scale the series for numerical stability
  scale_factor <- if (scaling && max(y) > 0) median(y[y > 0]) else 1
  y_scaled <- y / scale_factor
  occurrence <- as.numeric(y_scaled > 0)

  # Optimise parameters using Tweedie log-likelihood (p0 fixed, not estimated)
  opt <- twees_optimize(y_scaled, occurrence, damped)
  x <- opt$solution
  rho <- x[1]
  alpha_p <- x[2]
  theta_p <- if (damped) x[3] else 0
  mu0 <- x[4]
  alpha_mu <- x[5]
  theta_mu <- if (damped) x[6] else 0

  # Compute fitted values on the scaled series
  mu <- pmax(dampedSES(y_scaled, mu0, alpha_mu, theta_mu), .TWEES_EPSILON)
  p0 <- min(max(mean(occurrence), .TWEES_EPSILON), 1 - .TWEES_EPSILON)
  p <- pmin(pmax(dampedSES(occurrence, p0, alpha_p, theta_p), .TWEES_EPSILON), 1 - .TWEES_EPSILON)

  # Back-transform fitted values and residuals
  fitted <- mu * scale_factor
  residuals <- y - fitted

  structure(
    list(
      rho = rho,
      mu0 = mu0,
      alpha_mu = alpha_mu,
      theta_mu = theta_mu,
      p0 = p0,
      alpha_p = alpha_p,
      theta_p = theta_p,
      scale_factor = scale_factor,
      mean_y_scaled = mean(y_scaled),
      mean_occ = mean(occurrence),
      last_mu = mu[length(mu)],
      last_p = p[length(p)],
      last_y_scaled = y_scaled[length(y_scaled)],
      fitted = fitted,
      residuals = residuals
    ),
    class = "TWEES"
  )
}

#' Forecast a TWEES model
#'
#' Produces forecast distributions from a fitted TWEES
#' model using simulation.
#'
#' @inheritParams forecast.EMPDISTR
#' @param times The number of sample paths to use in estimating the forecast
#'   distribution.
#'
#' @return A distribution vector of forecasts: for h=1 the vector class is
#' `dist_tweedie`; for h>1 the vector class is `dist_sample`.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, TWEES(value))
#' forecast(fit, h = "7 days")
#'
#' @export
forecast.TWEES <- function(object, new_data, specials = NULL, times = 10000, ...) {
  h <- nrow(new_data)
  if (!is_integerish(times) || times <= 0) {
    abort("`times` must be a positive integer.")
  }

  # For the first step use a direct tweedie forecast
  mu_forecast <- object$alpha_mu * object$last_y_scaled +
    object$theta_mu * object$mean_y_scaled +
    (1 - object$alpha_mu - object$theta_mu) * object$last_mu
  mu_forecast <- max(mu_forecast, .TWEES_EPSILON)
  p_forecast <- object$alpha_p * as.integer(object$last_y_scaled > 0) +
    object$theta_p * object$mean_occ +
    (1 - object$alpha_p - object$theta_p) * object$last_p
  p_forecast <- min(max(p_forecast, .TWEES_EPSILON), 1 - .TWEES_EPSILON)
  phi <- -(mu_forecast^(2 - object$rho)) / ((2 - object$rho) * log1p(-p_forecast))
  dist_first <- dist_tweedie(
    mean = mu_forecast * object$scale_factor,
    dispersion = phi * object$scale_factor^(2 - object$rho),
    power = object$rho
  )

  if (h == 1) {
    return(dist_first)
  }
  sim <- twees_simulate(object, h, times)
  samples_rest <- as.list(as.data.frame(sim[, -1, drop = FALSE]))
  dist_rest <- dist_sample(samples_rest)

  c(dist_first, dist_rest)
}

#' Generate sample paths from a TWEES model
#'
#' @param x A fitted `TWEES` model object.
#' @inheritParams forecast.TWEES
#'
#' @return A vector of future paths from a dataset using a fitted model.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, TWEES(value))
#' generate(fit, new_data = tsibble::new_data(ts, 7))
#' @export
generate.TWEES <- function(x, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  sim <- twees_simulate(x, h, times = 1)
  new_data$.sim <- as.numeric(sim[1, ])
  new_data
}

#' Extract fitted values from a TWEES model
#'
#' @inherit fitted.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, TWEES(value))
#' fitted(fit)
#' @export
fitted.TWEES <- function(object, ...) {
  object$fitted
}

#' Extract residuals from a TWEES model
#'
#' @inherit residuals.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, TWEES(value))
#' residuals(fit)
#' @export
residuals.TWEES <- function(object, ...) {
  object$residuals
}


#' @export
model_sum.TWEES <- function(x) {
  if (x$theta_mu != 0) "TWEES(d)" else "TWEES(u)"
}

#' @export
tidy.TWEES <- function(x, ...) {
  terms <- c("alpha_mu", if (x$theta_mu != 0) "theta_mu",
             "alpha_p", if (x$theta_p != 0) "theta_p", "power", "mu0")
  ests  <- c(x$alpha_mu, if (x$theta_mu != 0) x$theta_mu,
             x$alpha_p, if (x$theta_p != 0) x$theta_p, x$rho, x$mu0)
  tibble(term = terms, estimate = ests)
}

#' @rdname TWEES
#' @export
report.TWEES <- function(object, ...) {
  cat("  Smoothing parameters:\n")
  cat(sprintf("    alpha_mu = %g\n", object$alpha_mu))
  if (object$theta_mu != 0) cat(sprintf("    theta_mu = %g\n", object$theta_mu))
  cat(sprintf("    alpha_p  = %g\n", object$alpha_p))
  if (object$theta_p != 0) cat(sprintf("    theta_p  = %g\n", object$theta_p))
  cat("\n  Tweedie parameters:\n")
  cat(sprintf("    power            = %g\n", object$rho))
  cat("\n  Initial state:\n")
  cat(sprintf("    mu[0] = %g\n", object$mu0))
  cat(sprintf("    p[0]  = %g  (fixed at mean(occurrence), not estimated)\n", object$p0))
  if (object$scale_factor != 1)
    cat(sprintf("\n  Scale factor: %g\n", object$scale_factor))
  invisible(object)
}

twees_simulate <- function(object, h, times) {
  forecast_samples <- matrix(NA_real_, nrow = times, ncol = h)

  # Build the state vector for the first-step mean (scaled)
  mu_state <- rep(
    object$alpha_mu * object$last_y_scaled +
      object$theta_mu * object$mean_y_scaled +
      (1 - object$alpha_mu - object$theta_mu) * object$last_mu,
    times
  )
  mu_state <- pmax(mu_state, .TWEES_EPSILON)
  p_state <- rep(
    object$alpha_p * as.integer(object$last_y_scaled > 0) +
      object$theta_p * object$mean_occ +
      (1 - object$alpha_p - object$theta_p) * object$last_p,
    times
  )
  p_state <- pmin(pmax(p_state, .TWEES_EPSILON), 1 - .TWEES_EPSILON)
  for (i in seq_len(h)) {
    phi <- -(mu_state^(2 - object$rho)) / ((2 - object$rho) * log1p(-p_state))
    y_new <- rtweedie(
      times,
      mean = mu_state,
      dispersion = phi,
      power = object$rho
    )
    forecast_samples[, i] <- y_new

    # Update the state on the scaled series
    mu_state <- object$alpha_mu * y_new +
      object$theta_mu * object$mean_y_scaled +
      (1 - object$alpha_mu - object$theta_mu) * mu_state
    mu_state <- pmax(mu_state, .TWEES_EPSILON)
    p_state <- object$alpha_p * as.integer(y_new > 0) +
      object$theta_p * object$mean_occ +
      (1 - object$alpha_p - object$theta_p) * p_state
    p_state <- pmin(pmax(p_state, .TWEES_EPSILON), 1 - .TWEES_EPSILON)
  }

  forecast_samples <- forecast_samples * object$scale_factor
  forecast_samples
}

twees_optimize <- function(y, occ, damped) {

  # Define the function to be optimised
  twees_nll <- function(x, y, occ) {
    rho <- x[1]
    alpha_p <- x[2]
    theta_p <- x[3]
    mu0 <- x[4]
    alpha_mu <- x[5]
    theta_mu <- x[6]

    mu <- dampedSES(y, mu0, alpha_mu, theta_mu)
    p0 <- min(max(mean(occ), .TWEES_EPSILON), 1 - .TWEES_EPSILON)
    p <- pmin(pmax(dampedSES(occ, p0, alpha_p, theta_p), .TWEES_EPSILON), 1 - .TWEES_EPSILON)
    phi <- -(mu^(2 - rho)) / ((2 - rho) * log1p(-p))
    -mean(dtweedie(y, mean = mu, dispersion = phi, power = rho, log = TRUE))
  }

  # In the undamped case set both damping parameters to 0
  if (!damped) {
    init_params <- c(1.5, 0.2, max(mean(y), .TWEES_EPSILON), 0.2)
    lb <- c(1.2 + .TWEES_EPSILON, rep(.TWEES_EPSILON, 3))
    ub <- c(1.8 - .TWEES_EPSILON, 1 - .TWEES_EPSILON, max(y) * 10, 1 - .TWEES_EPSILON)
    eval_f <- function(x) twees_nll(c(x[1], x[2], 0, x[3], x[4], 0), y, occ)
  } else {

    # Otherwise, learn both damping parameters with a simplex parametrisation
    init_params <- c(1.5, 0.2, 0.1 / (1 - 0.2), max(mean(y), .TWEES_EPSILON), 0.2, 0.1 / (1 - 0.2))
    lb <- c(1.2 + .TWEES_EPSILON, rep(.TWEES_EPSILON, 5))
    ub <- c(1.8 - .TWEES_EPSILON, rep(1 - .TWEES_EPSILON, 2), max(y) * 10, rep(1 - .TWEES_EPSILON, 2))
    eval_f <- function(x) twees_nll(c(x[1], x[2], (1 - x[2]) * x[3], x[4], x[5], (1 - x[5]) * x[6]), y, occ)
  }

  # Run the optimisation loop in an unconstrained way
  opt <- nloptr(
    x0 = init_params,
    eval_f = eval_f,
    lb = lb,
    ub = ub,
    opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
  )

  # Reparametrise the solution to the original parameterisation
  opt$solution <- if (damped) c(
    opt$solution[1], opt$solution[2], (1 - opt$solution[2]) * opt$solution[3],
    opt$solution[4:5], (1 - opt$solution[5]) * opt$solution[6]
  ) else c(opt$solution[1], opt$solution[2], 0, opt$solution[3:4], 0)
  opt
}

twees_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by TWEES.")
}
