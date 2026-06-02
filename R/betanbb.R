#' Beta-Negative Binomial Bayesian Dynamic Model
#'
#' Conjugate Bayesian dynamic model for count time series with a negative
#' binomial observation distribution and a Beta prior on the success probability
#' parameter. This extends the conjugate updating framework of Harvey &
#' Fernandes (1989) to the Beta-Negative Binomial family. The Beta prior is
#' updated at each time step using a discount factor `w`. Forecasts are available
#' as samples simulating from the model forward in time.
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
#' fc_ts <- ts |>
#'   model(BETANBB(value)) |>
#'   forecast(h = "7 days")
#'
#' fc_ts |> print()
#'
#'
#' if (requireNamespace("ggtime", quietly = TRUE)) {
#'   library(ggtime)
#'   fc_ts |> autoplot(ts)
#' }
#'
#' @importFrom fabletools new_model_class new_specials new_model_definition
#' @importFrom tsibble measured_vars
#' @importFrom rlang abort is_integerish
#' @importFrom distributional dist_sample
#' @importFrom nloptr nloptr
#' @importFrom stats rbeta rnbinom
#' @export
BETANBB <- function(formula, ...) {
  betanbb_model <- new_model_class(
    "BETANBB",
    train = train_betanbb,
    specials = new_specials(
      xreg = betanbb_no_xreg
    )
  )
  new_model_definition(betanbb_model, {{ formula }}, ...)
}

train_betanbb <- function(.data, specials, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by BETANBB.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by BETANBB.")
  }

  # Optimize parameters using negative log-likelihood
  opt <- betanbb_optimize(y)
  x <- opt$solution
  v <- x[1]
  a0 <- x[2]
  b0 <- x[3]
  w <- x[4]

  # Compute dynamic Beta parameters
  beta_params <- betaDynamic(y, v, a0, b0, w)
  a <- beta_params$a
  b <- beta_params$b

  # Compute expected fitted values with tower rule and residuals
  p_expected <- a / (a + b)
  fitted <- v * (1 - p_expected) / p_expected
  residuals <- y - fitted

  structure(
    list(
      v = v,
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
    class = "BETANBB"
  )
}

#' Forecast a BETANBB model
#'
#' Produces forecast distributions from a fitted BETANBB model using
#' simulation.
#'
#' @inheritParams forecast.EMPDISTR
#' @param times The number of sample paths to use in estimating the forecast
#'   distribution.
#'
#' @return A distribution vector of class `dist_sample`.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, BETANBB(value))
#' forecast(fit, h = "7 days")
#'
#' @importFrom fabletools forecast
#' @export
forecast.BETANBB <- function(object, new_data, specials = NULL, times = 10000, ...) {
  h <- nrow(new_data)
  if (!is_integerish(times) || times <= 0) {
    abort("`times` must be a positive integer.")
  }
  sim <- betanbb_simulate(object, h, times)
  samples <- as.list(as.data.frame(sim))
  dist_sample(samples)
}

#' Generate sample paths from a BETANBB model
#'
#' @param x A fitted `BETANBB` model object.
#' @inheritParams forecast.BETANBB
#'
#' @return A vector of future paths from a dataset using a fitted model.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, BETANBB(value))
#' generate(fit, new_data = tsibble::new_data(ts, 7))
#' @export
generate.BETANBB <- function(x, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  sim <- betanbb_simulate(x, h, 1L)
  new_data$.sim <- as.numeric(sim[1, ])
  new_data
}

#' Extract fitted values from a BETANBB model
#'
#' @inherit fitted.EMPDISTR
#'
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, BETANBB(value))
#' fitted(fit)
#' @export
fitted.BETANBB <- function(object, ...) {
  object$fitted
}

#' Extract residuals from a BETANBB model
#'
#' @inherit residuals.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, BETANBB(value))
#' residuals(fit)
#' @export
residuals.BETANBB <- function(object, ...) {
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
#' fit <- model(ts, BETANBB(value))
#' model_sum(fit[[1]][[1]])
#' @export
model_sum.BETANBB <- function(x) {
  "BETANBB"
}

betanbb_simulate <- function(object, h, times) {
  forecast_samples <- matrix(NA_real_, nrow = times, ncol = h)

  # Initialize Beta parameters with forward propagation of the last state
  a_state <- rep(
    object$w * object$last_a + (1 - object$w) + object$v,
    times
  )
  b_state <- rep(
    object$w * object$last_b + object$last_y,
    times
  )

  for (i in seq_len(h)) {
    # Sample p from Beta prior
    p_state <- rbeta(times, a_state, b_state)

    # Sample observations from Negative Binomial likelihood
    y_new <- rnbinom(times, object$v, p_state)
    forecast_samples[, i] <- y_new

    # Update Beta parameters
    a_state <- object$w * a_state + (1 - object$w) + object$v
    b_state <- object$w * b_state + y_new
  }

  forecast_samples
}


betanbb_optimize <- function(y) {

  # Define the negative log-likelihood function to be optimised
  nll_betanb <- function(x, y) {
    v <- x[1]
    a0 <- x[2]
    b0 <- x[3]
    w <- x[4]

    # Get the parameters following the updating rules
    beta_params <- betaDynamic(y, v, a0, b0, w)
    a <- beta_params$a
    b <- beta_params$b

    # Evaluate the negative log-likelihood of the model
    -mean(
      lbeta(v + a, y + b) - lbeta(v, y + 1) -
        lbeta(a, b) - log(v + y)
    )
  }

  # Run the optimization using nloptr with bounds
  nloptr(
    x0 = c(mean(y), 2, 2, 0.8),
    eval_f = function(x) nll_betanb(x, y),
    lb = c(.BETANBB_EPSILON, 1, .BETANBB_EPSILON, .BETANBB_EPSILON),
    ub = c(Inf, Inf, Inf, 1 - .BETANBB_EPSILON),
    opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
  )
}

betanbb_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by BETANBB.")
}

