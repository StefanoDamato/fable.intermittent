#' Gamma-Poisson Bayesian Dynamic Model
#'
#' Conjugate Bayesian dynamic model for count time series with a Poisson
#' observation distribution and a Gamma prior on the rate parameter, following
#' Harvey & Fernandes (1989). The Gamma prior is updated at each time step
#' using a discount factor `w` that controls how quickly past information
#' decays. The first-step forecast follows a Negative Binomial distribution,
#' and multi-step forecasts are obtained by simulating from the model forward in time.
#'
#' @param formula Model specification.
#' @param ... Not used.
#'
#' @references
#' Harvey, A. C., & Fernandes, C. (1989). Time series models for count or
#' qualitative observations. *Journal of Business & Economic Statistics*,
#' 7(4), 407--417. \doi{10.1080/07350015.1989.10509750}.
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
#'   model(GAMPOISB(value)) |>
#'   forecast(h = "7 days") |>
#'   ggtime::autoplot(ts)
#'
#' @importFrom fabletools new_model_class new_specials new_model_definition
#' @importFrom tsibble measured_vars
#' @importFrom rlang abort is_integerish
#' @importFrom distributional dist_sample dist_negative_binomial
#' @importFrom nloptr nloptr
#' @importFrom stats rgamma rpois
#' @export
GAMPOISB <- function(formula, ...) {
  gampoisb_model <- new_model_class(
    "GAMPOISB",
    train = train_gampoisb,
    specials = new_specials(
      xreg = gampoisb_no_xreg
    )
  )
  new_model_definition(gampoisb_model, {{ formula }}, ...)
}

train_gampoisb <- function(.data, specials, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by GAMPOISB.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by GAMPOISB.")
  }

  # Optimize parameters using negative log-likelihood
  opt <- gampoisb_optimize(y)
  x <- opt$solution
  a0 <- x[1]
  b0 <- x[2]
  w <- x[3]

  # Compute dynamic Gamma parameters
  gamma_params <- gammaDynamic(y, a0, b0, w)
  a <- gamma_params$a
  b <- gamma_params$b

  # Compute fitted values: E[Y | Gamma(a, b)] = a/b
  fitted <- a / b
  residuals <- y - fitted

  structure(
    list(
      a0 = a0,
      b0 = b0,
      w = w,
      a_state = a,
      b_state = b,
      last_y = y[length(y)],
      last_a = a[length(a)],
      last_b = b[length(b)],
      fitted = fitted,
      residuals = residuals
    ),
    class = "GAMPOISB"
  )
}

#' Forecast a GAMPOISB model
#'
#' Produces forecast distributions from a fitted GAMPOISB model using
#' simulation.
#'
#' @inheritParams forecast.EMPDISTR
#' @param times The number of sample paths to use in estimating the forecast
#'   distribution.
#'
#' @return A distribution vector of forecasts: for h=1 the vector class is
#' `dist_negative_binomial`; for h>1 the vector class is `dist_sample`.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, GAMPOISB(value))
#' forecast(fit, h = "7 days")
#'
#' @export
forecast.GAMPOISB <- function(object, new_data, specials = NULL, times = 10000, ...) {
  h <- nrow(new_data)
  if (!is_integerish(times) || times <= 0) {
    abort("`times` must be a positive integer.")
  }

  # Initialize Gamma parameters with forward propagation
  a_forecast <- object$w * object$last_a + object$last_y
  b_forecast <- object$w * object$last_b + 1
  dist_first <- dist_negative_binomial(size = a_forecast, prob = b_forecast / (b_forecast + 1))

  if (h == 1) {
    return(dist_first)
  }

  sim <- gampoisb_simulate(object, h, times)
  samples_rest <- as.list(as.data.frame(sim[, -1, drop = FALSE]))
  dist_rest <- dist_sample(samples_rest)

  c(dist_first, dist_rest)
}

#' Generate sample paths from a GAMPOISB model
#'
#' @param x A fitted `GAMPOISB` model object.
#' @inheritParams forecast.GAMPOISB
#'
#' @return A vector of future paths from a dataset using a fitted model.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, GAMPOISB(value))
#' generate(fit, new_data = tsibble::new_data(ts, 7))
#' @export
generate.GAMPOISB <- function(x, new_data, specials = NULL, ...) {
  h <- NROW(new_data)
  sim <- gampoisb_simulate(x, h, 1L)
  new_data$.sim <- as.numeric(sim[1, ])
  new_data
}

#' Extract fitted values from a GAMPOISB model
#'
#' @inherit fitted.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, GAMPOISB(value))
#' fitted(fit)
#' @export
fitted.GAMPOISB <- function(object, ...) {
  object$fitted
}

#' Extract residuals from a GAMPOISB model
#'
#' @inherit residuals.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, GAMPOISB(value))
#' residuals(fit)
#' @export
residuals.GAMPOISB <- function(object, ...) {
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
#' fit <- model(ts, GAMPOISB(value))
#' model_sum(fit[[1]][[1]])
#' @export
model_sum.GAMPOISB <- function(x) {
  "GAMPOISB"
}

gampoisb_simulate <- function(object, h, times) {
  forecast_samples <- matrix(NA_real_, nrow = times, ncol = h)

  # Initialize Gamma parameters with forward propagation
  a_state <- rep(
    object$w * object$last_a + object$last_y,
    times
  )
  b_state <- rep(
    object$w * object$last_b + 1,
    times
  )

  for (i in seq_len(h)) {
    # Sample lambda from Gamma prior
    lambda_state <- rgamma(times, a_state, b_state)

    # Sample observations from Poisson likelihood
    y_new <- rpois(times, lambda_state)
    forecast_samples[, i] <- y_new

    # Update Gamma parameters
    a_state <- object$w * a_state + y_new
    b_state <- object$w * b_state + 1
  }

  forecast_samples
}

gampoisb_optimize <- function(y) {

  # Define the negative log-likelihood function to be optimised
  nll_gampois <- function(x, y) {
    a0 <- x[1]
    b0 <- x[2]
    w <- x[3]

    # Compute the dynamic parameters
    gamma_params <- gammaDynamic(y, a0, b0, w)
    a <- gamma_params$a
    b <- gamma_params$b

    # Evaluate the negative log-likelihood of the model
    -mean(
      lchoose(a + y - 1, y) +
        a * log(b) - (a + y) * log(1 + b)
    )
  }

  # Run the optimization using nloptr with bounds
  nloptr(
    x0 = c(1, 1, 0.8),
    eval_f = function(x) nll_gampois(x, y),
    lb = c(.GAMPOISB_EPSILON, .GAMPOISB_EPSILON, .GAMPOISB_EPSILON),
    ub = c(Inf, Inf, 1 - .GAMPOISB_EPSILON),
    opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
  )
}

gampoisb_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by GAMPOISB.")
}

