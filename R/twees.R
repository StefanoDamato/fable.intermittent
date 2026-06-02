#' Tweedie Exponential Smoothing
#'
#' Exponential smoothing state space model for intermittent demand with a
#' Tweedie observation distribution. The conditional mean of the Tweedie
#' is governed by a (optionally damped) exponential smoothing process.
#' The power and dispersion parameter are estimated to maximise the likelihood.
#' The Tweedie family naturally models both zeros and large spikes
#' via its compound Poisson-Gamma nature. The first-step forecast follows
#' a Tweedie distribution, and multi-step forecasts are obtained by simulating
#' from the model forward in time. The model parameters are estimated by
#'
#' @param formula Model specification.
#' @param damped Logical. If `TRUE` (default), the exponential smoothing
#'   component uses a damping parameter.
#' @param scaling Logical. If `TRUE` (default), the time series is divided by
#'   its maximum value before fitting and predictions are back-transformed.
#'   This improves numerical stability.
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

  # Optimise parameters using Tweedie log-likelihood
  opt <- twees_optimize(y_scaled, damped)
  x <- opt$solution
  phi <- x[1]
  power <- x[2]
  mu0 <- x[3]
  alpha <- x[4]
  theta <- if (damped) x[5] else 0

  # Compute fitted values on the scaled series
  mu <- dampedSES(y_scaled, mu0, alpha, theta)
  mu <- pmax(mu, .TWEES_EPSILON)

  # Back-transform fitted values and residuals
  fitted <- mu * scale_factor
  residuals <- y - fitted

  structure(
    list(
      phi = phi,
      power = power,
      mu0 = mu0,
      alpha = alpha,
      theta = theta,
      scale_factor = scale_factor,
      mean_y_scaled = mean(y_scaled),
      last_mu = mu[length(mu)],
      last_y_scaled = y_scaled[length(y_scaled)],
      fitted = fitted,
      residuals = residuals
    ),
    class = "TWEES"
  )
}

#' Forecast a TWEES model
#'
#' Produces forecast distributions from a fitted TWEES model using simulation.
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
  mu_forecast <- object$alpha * object$last_y_scaled +
    object$theta * object$mean_y_scaled +
    (1 - object$alpha - object$theta) * object$last_mu
  mu_forecast <- max(mu_forecast, .TWEES_EPSILON)
  dist_first <- dist_tweedie(
    mean = mu_forecast * object$scale_factor,
    dispersion = object$phi * object$scale_factor^(2 - object$power),
    power = object$power
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


#' @inherit model_sum.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, TWEES(value))
#' model_sum(fit[[1]][[1]])
#' @export
model_sum.TWEES <- function(x) {
  "TWEES"
}

twees_simulate <- function(object, h, times) {
  forecast_samples <- matrix(NA_real_, nrow = times, ncol = h)

  # Build the state vector for the first-step mean (scaled)
  mu_state <- rep(
    object$alpha * object$last_y_scaled +
      object$theta * object$mean_y_scaled +
      (1 - object$alpha - object$theta) * object$last_mu,
    times
  )
  mu_state <- pmax(mu_state, .TWEES_EPSILON)

  for (i in seq_len(h)) {
    # Sample from Tweedie on the original scale
    y_new <- rtweedie(
      times,
      mean = mu_state,
      dispersion = object$phi,
      power = object$power
    )
    forecast_samples[, i] <- y_new

    # Update the state on the scaled series
    mu_state <- object$alpha * y_new +
      object$theta * object$mean_y_scaled +
      (1 - object$alpha - object$theta) * mu_state
    mu_state <- pmax(mu_state, .TWEES_EPSILON)
  }

  forecast_samples <- forecast_samples * object$scale_factor
  forecast_samples
}

twees_optimize <- function(y, damped) {

  # Define the function to be optimised
  twees_nll <- function(x, y) {
    phi <- x[1]
    power <- x[2]
    mu0 <- x[3]
    alpha <- x[4]
    theta <- x[5]

    # Fit the exponential smoothing and return the negative log-likkelihood
    mu <- dampedSES(y, mu0, alpha, theta)
    -mean(dtweedie(y, mean = mu, dispersion = phi, power = power, log = TRUE))
  }

  # Define good starting values based on moments
  mean_y <- mean(y)
  var_y <- var(y)
  max_y <- max(y[y > 0], na.rm = TRUE)
  rho_init <- 1.5
  phi_init <- (mean_y^rho_init) / var_y
  phi_upper <- 10 * max(mean_y, mean_y^2) / var_y

  # In the undamped case specify the parameter vector with theta fixed to 0
  if (!damped) {
    init_params <- c(phi_init, rho_init, max(mean_y, .TWEES_EPSILON), 0.3)
    lb <- c(.TWEES_EPSILON, 1 + .TWEES_EPSILON, .TWEES_EPSILON, .TWEES_EPSILON)
    ub <- c(phi_upper, 2 - .TWEES_EPSILON, max_y * 10, 1 - .TWEES_EPSILON)

    # Run the optimistion with bounds using nloptr
    opt <- nloptr(
      x0 = init_params,
      eval_f = function(x) twees_nll(c(x, 0), y),
      lb = lb,
      ub = ub,
      opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
    )
    opt$solution <- c(opt$solution, 0)
  } else {

    # In the damped case, specify the full parameter vector
    init_params <- c(phi_init, rho_init, max(mean_y, .TWEES_EPSILON), 0.3, 0.1)
    lb <- c(.TWEES_EPSILON, 1 + .TWEES_EPSILON, .TWEES_EPSILON, .TWEES_EPSILON, 0)
    ub <- c(phi_upper, 2 - .TWEES_EPSILON, max_y * 10, 1 - .TWEES_EPSILON, 1)

    # Run the optimization with bounds and a linear constraint using nloptr
    opt <- nloptr(
      x0 = init_params,
      eval_f = function(x) twees_nll(x, y),
      lb = lb,
      ub = ub,
      eval_g_ineq = function(x) x[4] + x[5] - 1 + .TWEES_EPSILON,
      opts = list(algorithm = "NLOPT_LN_COBYLA", maxeval = 500)
    )
  }

  opt
}

twees_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by TWEES.")
}

