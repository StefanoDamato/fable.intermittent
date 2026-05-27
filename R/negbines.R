#' Negative Binomial Exponential Smoothing
#'
#' Exponential smoothing state space model for intermittent demand with a
#' negative binomial observation distribution, as proposed by Snyder, Ord &
#' Beaumont (2012). The conditional mean of the negative binomial is governed by
#' a (optionally damped) exponential smoothing process. The probability
#' parameter is estimated to maximise the likelihood. The first-step forecast 
#' follows a Negative Binomial distribution, and multi-step forecasts 
#' are obtained by simulating from the model forward in time.
#'
#' @param formula Model specification.
#' @param damped Logical. If `TRUE` (default), the exponential smoothing
#'   component uses a damping parameter.
#' @param ... Not used.
#'
#' @references
#' Snyder, R. D., Ord, J. K., & Beaumont, A. (2012). Forecasting the
#' intermittent demand for slow-moving inventories: A modelling approach.
#' *International Journal of Forecasting*, 28(2), 485--496.
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
#'   model(NEGBINES(value)) |>
#'   forecast(h = "7 days") |>
#'   ggtime::autoplot(ts)
#'
#' @importFrom fabletools new_model_class new_specials new_model_definition
#' @importFrom tsibble measured_vars
#' @importFrom rlang abort is_integerish
#' @importFrom distributional dist_sample dist_negative_binomial
#' @importFrom nloptr nloptr
#' @importFrom stats rnbinom dnbinom
#' @export
NEGBINES <- function(formula, damped = TRUE, ...) {
  negbines_model <- new_model_class(
    "NEGBINES",
    train = train_negbines,
    specials = new_specials(
      xreg = negbines_no_xreg
    )
  )
  new_model_definition(negbines_model, {{ formula }}, damped = damped, ...)
}

train_negbines <- function(.data, specials, damped, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by NEGBINES.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by NEGBINES.")
  }

  if (!is.logical(damped)) {
    abort("`damped` must be a boolean.")
  }

  # Optimize parameters using Negative Binomial likelihood 
  opt <- negbines_optimize(y, damped)
  x <- opt$solution
  prob <- x[1]
  mu0 <- x[2]
  alpha <- x[3]
  phi <- ifelse(damped, x[4], 0)

  # Compute fitted values and residuals
  mu <- dampedSES(y, mu0, alpha, phi)
  fitted <- mu
  residuals <- y - fitted

  structure(
    list(
      prob = prob,
      mu0 = mu0,
      alpha = alpha,
      phi = phi,
      mean_y = mean(y),
      last_mu = mu[length(mu)],
      last_y = y[length(y)],
      fitted = fitted,
      residuals = residuals
    ),
    class = "NEGBINES"
  )
}

#' Forecast a NEGBINES model
#'
#' Produces forecast distributions from a fitted NEGBINES model using
#' simulation.
#'
#' @inheritParams forecast.EMPDISTR
#' @param times The number of sample paths to use in estimating the forecast
#'   distribution.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, NEGBINES(value))
#' forecast(fit, h = "7 days")
#'
#' @export
forecast.NEGBINES <- function(object, new_data, specials = NULL, times = 10000, ...) {
  h <- nrow(new_data)
  if (!is_integerish(times) || times <= 0) {
    abort("`times` must be a positive integer.")
  }

  # First step is always direct Negative Binomial
  mu_forecast <- object$alpha * object$last_y +
    object$phi * object$mean_y +
    (1 - object$alpha - object$phi) * object$last_mu
  
  size <- mu_forecast * object$prob / (1 - object$prob)
  dist_first <- dist_negative_binomial(size = size, prob = object$prob)
  
  if (h == 1) {
    return(dist_first)
  }
  
  sim <- negbines_simulate(object, h, times)
  samples_rest <- as.list(as.data.frame(sim[, -1, drop = FALSE]))
  dist_rest <- dist_sample(samples_rest)
  
  c(dist_first, dist_rest)
}

#' Extract fitted values from a NEGBINES model
#'
#' @inheritParams forecast.NEGBINES
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, NEGBINES(value))
#' fitted(fit)
#' @export
fitted.NEGBINES <- function(object, ...) {
  object$fitted
}

#' Extract residuals from a NEGBINES model
#'
#' @inheritParams forecast.NEGBINES
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, NEGBINES(value))
#' residuals(fit)
#' @export
residuals.NEGBINES <- function(object, ...) {
  object$residuals
}

#' @export
model_sum.NEGBINES <- function(x) {
  "NEGBINES"
}

#' Generate sample paths from a NEGBINES model
#'
#' @param x A fitted `NEGBINES` model object.
#' @inheritParams forecast.NEGBINES
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, NEGBINES(value))
#' generate(fit, new_data = tsibble::new_data(ts, 7))
#' @export
generate.NEGBINES <- function(x, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  sim <- negbines_simulate(x, h, times = 1)
  new_data$.sim <- as.numeric(sim[1, ])
  new_data
}

negbines_simulate <- function(object, h, times) {
  forecast_samples <- matrix(NA_real_, nrow = times, ncol = h)

  # Build the state vector for the firststep mean
  mu_state <- rep(
    object$alpha * object$last_y +
      object$phi * object$mean_y +
      (1 - object$alpha - object$phi) * object$last_mu,
    times
  )

  # Sequentially update the processes and sample forecasts
  for (i in seq_len(h)) {
    y_new <- rnbinom(
      times,
      size = mu_state * object$prob / (1 - object$prob),
      prob = object$prob
    )
    forecast_samples[, i] <- y_new

    # Update the state vector for the next step
    mu_state <- object$alpha * y_new +
      object$phi * object$mean_y +
      (1 - object$alpha - object$phi) * mu_state
  }

  forecast_samples
}

negbines_optimize <- function(y, damped) {

  # Define the negative log-likelihood function to be optimised
  negbines_nll <- function(x, y) {
    prob <- x[1]
    mu0 <- x[2]
    alpha <- x[3]
    phi <- x[4]

    # Fit the damped exponential smoothing and return the negative log-likelihood
    mu <- dampedSES(y, mu0, alpha, phi)
    -mean(dnbinom(
      y,
      size = mu * prob / (1 - prob),
      prob = prob,
      log = TRUE
    ))
  }
  
  # In the undamped case set the last parameter to 0
  if (!damped) {
    init_params <- c(0.5, mean(y), 0.3)
    lb <- c(negbines_epsilon, negbines_epsilon, negbines_epsilon)
    ub <- c(1 - negbines_epsilon, max(y) * 10, 1 - negbines_epsilon)
    
    # Run the bounded optimization using nloptr
    opt <- nloptr(
      x0 = init_params,
      eval_f = function(x) negbines_nll(c(x, 0), y),
      lb = lb,
      ub = ub,
      opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
    )
    opt$solution <- c(opt$solution, 0)
  } else {

    # In the damped case specify the full parameter vector
    init_params <- c(0.5, mean(y), 0.3, 0.1)
    lb <- c(negbines_epsilon, negbines_epsilon, negbines_epsilon, 0)
    ub <- c(1 - negbines_epsilon, max(y), 1 - negbines_epsilon, 1)
    
    # run the bounded optimization with a linear constraint using nloptr
    opt <- nloptr(
      x0 = init_params,
      eval_f = function(x) negbines_nll(x, y),
      lb = lb,
      ub = ub,
      eval_g_ineq = function(x) x[3] + x[4] - 1 + negbines_epsilon,
      opts = list(algorithm = "NLOPT_LN_COBYLA", maxeval = 500)
    )
  }
  
  opt
}


negbines_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by NEGBINES.")
}

negbines_epsilon <- 1e-4
