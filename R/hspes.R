#' Hurdle-Shifted Poisson Exponential Smoothing
#'
#' Exponential smoothing state space model for intermittent demand with a
#' hurdle-shifted Poisson observation distribution, following the framework of
#' Snyder, Ord & Beaumont (2012). The model decomposes demand into two
#' components: a Bernoulli occurrence process (probability of non-zero demand)
#' and a shifted Poisson demand-size process. Both components are driven by
#' independent (optionally damped) exponential smoothing state equations.
#' The first-step forecast follows an hurdle-shifted Poisson distribution,
#' multi-step forecasts are obtained by simulating from the model forward in time.
#'
#' @param formula Model specification.
#' @param damped Logical. If `TRUE` (default), both the occurrence and demand
#'   smoothing components use a damping parameter.
#' @param ... Not used.
#'
#' @references
#' Snyder, R. D., Ord, J. K., & Beaumont, A. (2012). Forecasting the
#' intermittent demand for slow-moving inventories: A modelling approach.
#' *International Journal of Forecasting*, 28(2), 485--496.
#' \doi{10.1016/j.ijforecast.2011.03.009}.
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
#' ts |>
#'   model(HSPES(value)) |>
#'   forecast(h = "7 days") |>
#'   ggtime::autoplot(ts)
#'
#' @importFrom fabletools new_model_class new_specials new_model_definition
#' @importFrom tsibble measured_vars
#' @importFrom rlang abort is_integerish
#' @importFrom distributional dist_sample dist_poisson dist_transformed dist_inflated
#' @importFrom nloptr nloptr
#' @importFrom stats runif rpois dbinom dpois
#' @importFrom utils tail
#' @export
HSPES <- function(formula, damped = TRUE, ...) {
  hspes_model <- new_model_class(
    "HSPES",
    train = train_hspes,
    specials = new_specials(
      xreg = hspes_no_xreg
    )
  )
  new_model_definition(hspes_model, {{ formula }}, damped = damped, ...)
}

#' @importFrom utils tail
train_hspes <- function(.data, specials, damped, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by HSPES.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by HSPES.")
  }
  if (all(y == 0)) {
    abort("The time series is all zero.")
  }

  if (!is.logical(damped)) {
    abort("`damped` must be a boolean.")
  }

  # Croston's decomposition
  decomp <- crostons_decomp(y)
  occurrence <- decomp$occurrence
  intervals <- decomp$intervals
  shifted_demand <- decomp$demand - 1

  # Optimize occurrence and demand components separately
  opt_occ <- hspes_optimize_occurrence(occurrence, damped)
  opt_dem <- hspes_optimize_demand(shifted_demand, damped)

  # Extract parameters for occurrence
  x_occ <- opt_occ$solution
  p0 <- x_occ[1]
  alpha_occ <- x_occ[2]
  phi_occ <- ifelse(damped, x_occ[3], 0)

  # Extract parameters for demand
  x_dem <- opt_dem$solution
  lambda0 <- x_dem[1]
  alpha_dem <- x_dem[2]
  phi_dem <- ifelse(damped, x_dem[3], 0)

  # Compute fitted exponential smoothing states
  p <- dampedSES(occurrence, p0, alpha_occ, phi_occ)
  lambda <- if (length(shifted_demand) == 1) {
    rep(lambda0, length(shifted_demand))
  } else {
    dampedSES(shifted_demand, lambda0, alpha_dem, phi_dem)
  }

  # Compute fitted values and residuals
  tail_len <- length(y) - sum(intervals)
  lambda_pred <- alpha_dem * tail(shifted_demand, 1) +
    phi_dem * mean(shifted_demand) +
    (1 - alpha_dem - phi_dem) * tail(lambda, 1)
  lambda_fitted <- rep(
    c(lambda, lambda_pred),
    times = c(intervals, tail_len)
  )
  fitted <- p * (lambda_fitted + 1)
  residuals <- y - fitted

  # Save model components in a structured object
  structure(
    list(
      p0 = p0,
      alpha_occ = alpha_occ,
      phi_occ = phi_occ,
      lambda0 = lambda0,
      alpha_dem = alpha_dem,
      phi_dem = phi_dem,
      mean_occurrence = mean(occurrence),
      mean_shifted_demand = mean(shifted_demand),
      last_occurrence = occurrence[length(occurrence)],
      last_shifted_demand = shifted_demand[length(shifted_demand)],
      last_p = p[length(p)],
      last_lambda = lambda[length(lambda)],
      fitted = fitted,
      residuals = residuals
    ),
    class = "HSPES"
  )
}

#' Forecast a HSPES model
#'
#' Produces forecast distributions from a fitted HSPES model using
#' simulation.
#'
#' @inheritParams forecast.EMPDISTR
#' @param times The number of sample paths to use in estimating the forecast
#'   distribution.
#'
#' @return A distribution vector of forecasts: for h=1 the vector class is
#' `dist_inflated` (hurdle-shifted Poisson); for h>1 the vector class is `dist_sample`.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, HSPES(value))
#' forecast(fit, h = "7 days")
#'
#' @export
forecast.HSPES <- function(object, new_data, specials = NULL, times = 10000, ...) {
  h <- nrow(new_data)
  if (!is_integerish(times) || times <= 0) {
    abort("`times` must be a positive integer.")
  }

  # Predict the next p and lambda values
  p_forecast <- object$alpha_occ * object$last_occurrence +
    object$phi_occ * object$mean_occurrence +
    (1 - object$alpha_occ - object$phi_occ) * object$last_p

  # Predict the next lambda value
  lambda_forecast <- object$alpha_dem * object$last_shifted_demand +
    object$phi_dem * object$mean_shifted_demand +
    (1 - object$alpha_dem - object$phi_dem) * object$last_lambda

  # Build the distributional hurdle-shifted Poisson for h = 1
  dist_hsp <- dist_poisson(lambda_forecast) |>
    dist_transformed(\(x) x + 1, \(x) x - 1) |>
    dist_inflated(1 - p_forecast)
  if (h == 1) {
    return(dist_hsp)
  }

  # For h > 1 include sampled distributions from step 2 onwards
  sim <- hspes_simulate(object, h, times)
  samples_rest <- as.list(as.data.frame(sim[, -1, drop = FALSE]))
  dist_rest <- dist_sample(samples_rest)

  c(dist_hsp, dist_rest)
}

#' Generate sample paths from a HSPES model
#'
#' @param x A fitted `HSPES` model object.
#' @inheritParams forecast.HSPES
#'
#' @return A vector of future paths from a dataset using a fitted model.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, HSPES(value))
#' generate(fit, new_data = tsibble::new_data(ts, 7))
#' @export
generate.HSPES <- function(x, new_data, specials = NULL, ...) {
  h <- NROW(new_data)
  sim <- hspes_simulate(x, h, 1L)
  new_data$.sim <- as.numeric(sim[1, ])
  new_data
}

#' Extract fitted values from a HSPES model
#'
#' @inherit fitted.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, HSPES(value))
#' fitted(fit)
#' @export
fitted.HSPES <- function(object, ...) {
  object$fitted
}

#' Extract residuals from a HSPES model
#'
#' @inherit residuals.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, HSPES(value))
#' residuals(fit)
#' @export
residuals.HSPES <- function(object, ...) {
  object$residuals
}

#' Return model name
#' @inherit model_sum.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, HSPES(value))
#' model_sum(fit[[1]][[1]])
#' @export
model_sum.HSPES <- function(x) {
  "HSPES"
}

hspes_simulate <- function(object, h, times) {
  forecast_samples <- matrix(0, nrow = times, ncol = h)

  # Build the state vector p
  p_state <- rep(
    object$alpha_occ * object$last_occurrence +
      object$phi_occ * object$mean_occurrence +
      (1 - object$alpha_occ - object$phi_occ) * object$last_p,
    times
  )

  # Build the state vector for lambda
  lambda_state <- rep(
    object$alpha_dem * object$last_shifted_demand +
      object$phi_dem * object$mean_shifted_demand +
      (1 - object$alpha_dem - object$phi_dem) * object$last_lambda,
    times
  )

  # Sequentially update the processes and sample forecasts
  for (i in seq_len(h)) {
    occ_new <- as.integer(runif(times) <= p_state)
    dem_new <- rpois(times, lambda_state)
    forecast_samples[, i] <- ifelse(occ_new == 1, dem_new + 1, 0)

    # Update states
    p_state <- object$alpha_occ * occ_new +
      object$phi_occ * object$mean_occurrence +
      (1 - object$alpha_occ - object$phi_occ) * p_state
    lambda_state <- ifelse(occ_new == 1,
      object$alpha_dem * dem_new +
        object$phi_dem * object$mean_shifted_demand +
        (1 - object$alpha_dem - object$phi_dem) * lambda_state,
      lambda_state
    )
  }

  forecast_samples
}

hspes_optimize_occurrence <- function(occurrence, damped) {

  # Define the bernoulli negative log-likelihood function to be optimised
  nll_occ <- function(x, occurrence) {
    p0 <- x[1]
    alpha_occ <- x[2]
    phi_occ <- x[3]

    # Fit the damped exponential smoothing and return the negative
    p <- dampedSES(occurrence, p0, alpha_occ, phi_occ)
    -mean(dbinom(occurrence, size = 1, prob = p, log = TRUE))
  }

  # In the undamped case, set the last parameter to 0
  if (!damped) {
    init_params <- c(min(mean(occurrence), 1 - .HSPES_EPSILON), 0.2)
    lb <- c(.HSPES_EPSILON, .HSPES_EPSILON)
    ub <- c(1 - .HSPES_EPSILON, 1 - .HSPES_EPSILON)

    # Run the optimization using nloptr with bounds
    opt <- nloptr(
      x0 = init_params,
      eval_f = function(x) nll_occ(c(x, 0), occurrence),
      lb = lb,
      ub = ub,
      opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
    )
  } else {

    # In the damped case, specify the full parameter vector
    init_params <- c(min(mean(occurrence), 1 - .HSPES_EPSILON), 0.2, 0.2)
    lb <- c(.HSPES_EPSILON, .HSPES_EPSILON, .HSPES_EPSILON)
    ub <- c(1 - .HSPES_EPSILON, 1 - .HSPES_EPSILON, 1 - .HSPES_EPSILON)

    # Run the optimization using nloptr with bounds and a linear constraint
    opt <-nloptr(
      x0 = init_params,
      eval_f = function(x) nll_occ(x, occurrence),
      lb = lb,
      ub = ub,
      eval_g_ineq = function(x) x[2] + x[3] - 1 + .HSPES_EPSILON,
      opts = list(algorithm = "NLOPT_LN_COBYLA", maxeval = 500)
    )
  }

  opt
}

hspes_optimize_demand <- function(shifted_demand, damped) {

  # Define the poisson negative log-likelihood function to be optimised
  nll_dem <- function(x, shifted_demand) {
    lambda0 <- x[1]
    alpha_dem <- x[2]
    phi_dem <- x[3]

    # Fit the damped exponential smoothing and return the negative log-likelihood
    lambda <- dampedSES(shifted_demand, lambda0, alpha_dem, phi_dem)
    -mean(dpois(shifted_demand, lambda, log = TRUE))
  }

  # When demand is constant optimisation is not needed
  if (length(unique(shifted_demand)) == 1) {
    lambda0 <- max(unique(shifted_demand), .HSPES_EPSILON)
    opt <- list(solution = c(lambda0, 0, 0))
    return(opt)
  }

  # In the undamped case, set the last parameter to 0
  if (!damped) {
    init_params <- c(max(mean(shifted_demand), .HSPES_EPSILON), 0.2)
    lb <- c(.HSPES_EPSILON, .HSPES_EPSILON)
    ub <- c(max(shifted_demand) * 10, 1 - .HSPES_EPSILON)

    # Run the optimization using nloptr with bounds
    opt <-nloptr(
      x0 = init_params,
      eval_f = function(x) nll_dem(c(x, 0), shifted_demand),
      lb = lb,
      ub = ub,
      opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
    )
    opt$solution <- c(opt$solution, 0)
  } else {

    # In the damped case, specify the full parameter vector
    init_params <- c(max(mean(shifted_demand), .HSPES_EPSILON), 0.2, 0.2)
    lb <- c(.HSPES_EPSILON, .HSPES_EPSILON, .HSPES_EPSILON)
    ub <- c(max(shifted_demand) * 10, 1 - .HSPES_EPSILON, 1 - .HSPES_EPSILON)

    # run the optimization using nloptr with bounds and a linear constraint
    opt <- nloptr(
      x0 = init_params,
      eval_f = function(x) nll_dem(x, shifted_demand),
      lb = lb,
      ub = ub,
      eval_g_ineq = function(x) x[2] + x[3] - 1 + .HSPES_EPSILON,
      opts = list(algorithm = "NLOPT_LN_COBYLA", maxeval = 500)
    )
  }

  opt
}

hspes_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by HSPES.")
}

