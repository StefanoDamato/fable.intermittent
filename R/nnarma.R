#' Non-negative ARMA model
#'
#' A non-negative autoregressive moving average (NNARMA) model for
#' intermittent demand, where the demand is driven by an ARMA(1, 1), constrained
#' to remain non-negative. Forecasts are returned as Gaussian distributions 
#' for each time step, with the negative probability mass collapsed to zero.
#'
#' @param formula Model specification.
#' @param object A fitted model object.
#' @param ... Not used.
#'
#' @references
#'
#' Sbrana, G., & Babai, M. Z. (2026). Non-negative autoregressive moving average
#' models for intermittent demand: Forecast accuracy and inventory implications.
#' *European Journal of Operational Research*. \doi{10.1016/j.ejor.2026.06.009}.
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
#'   model(NNARMA(value)) |>
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
#' @export
NNARMA <- function(formula, ...) {
  nnarma_model <- new_model_class(
    "NNARMA",
    train = train_nnarma,
    specials = new_specials(
      xreg = nnarma_no_xreg
    )
  )
  new_model_definition(nnarma_model, {{ formula }}, ...)
}

train_nnarma <-function(.data, specials, ...) {
  y <- unclass(.data)[[measured_vars(.data)]]
  
  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by NNARMA.")
  }
  if (all(y == 0)) {
    abort("The time series is all zero.")
  }
  
  # Extract the correct frequency and deseasonalize the data
  max_prop_zeros <- 0.95
  period <- get_freq(.data)
  if (period < 1) {
    abort("The seasonal period must be greater than or equal to 1.")
  }
  deseasonalized <- deseasonalize(y, period = period,
                                  max_prop_zeros = max_prop_zeros,
                                  normalize = TRUE)
  y_deseasonalized <- deseasonalized$y_deseasonalized
  seasons <- deseasonalized$seasons
  
  opt <- nnarma_optimize(y_deseasonalized)
  x <- opt$solution
  phi <- x[1]
  theta <- x[2]
  co <- x[3]
  
  process <- armaDynamic(y_deseasonalized, phi, theta, co)
  v <- process$v
  m <- process$m

  # Deseasonalize fitted values and residuals
  fitted <- (m + co) * if (is.null(seasons)) 1 else rep(seasons, length.out = length(y))
  residuals <- y - fitted
  
  structure(
    list(
      phi = phi,
      theta = theta,
      co = co,
      frequency = period,
      seasons = seasons,
      v_state = v,
      last_v = v[length(v)],
      last_m = m[length(m)],
      last_y = y[length(y)],
      fitted = fitted,
      residuals = residuals
    ),
    class = "NNARMA"
  )
}

#' Forecast a NNARMA model
#'
#' Produces forecast distributions from a fitted NNARMA model.
#'
#' @inheritParams forecast.EMPSD
#'
#' @return A distribution vector of class `dist_normal_nonneg`.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, NNARMA(value))
#' forecast(fit, h = "7 days")
#'
#' @export
forecast.NNARMA <- function(object, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  
  # Compute the forecast mean and variance 
  K <- object$phi + object$theta
  var_v <- var(object$v_state)
  if (!is.finite(var_v) || var_v < .NNARMA_EPSILON) {
    var_v <- .NNARMA_EPSILON
  }
  pows <- 2 * (seq_len(h - 1) - 1)
  mean_fc <- (object$phi * object$last_m + K * object$last_v) * object$phi^(0:(h - 1)) +
    object$co
  var_fc <- c(0, K^2 * var_v * cumsum(object$phi^pows)) + var_v

  # Adjust for seasonality if necessary
  if (!is.null(object$seasons)) {
    s <- 1 + length(object$v_state) %% object$frequency
    seasons <- rep(object$seasons,
                   1 + ceiling((h - object$frequency + s - 1) / object$frequency))
    mean_fc <- mean_fc * seasons[1:h + s -1]
  }
  
  # Return the Gaussian forecast distribution (floor sd to avoid degenerate N(0,0))
  sd_fc <- pmax(sqrt(pmax(var_fc, .NNARMA_EPSILON)), .NNARMA_EPSILON)
  dist_normal_nonneg(mean_fc, sd_fc)
}

#' Extract fitted values from a NNARMA model
#'
#' @inherit fitted.EMPSD
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, NNARMA(value))
#' fitted(fit)
#'
#' @export
fitted.NNARMA <- function(object, ...) {
  object$fitted
}

#' Extract residuals from a NNARMA model
#'
#' @inherit residuals.EMPSD
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, NNARMA(value))
#' residuals(fit)
#' @export
residuals.NNARMA <- function(object, ...) {
  object$residuals
}


#' @export
model_sum.NNARMA <- function(x) {
  "NNARMA"
}

#' @importFrom generics tidy
#' @importFrom tibble tibble
#' @export
tidy.NNARMA <- function(x, ...) {
  tibble(
    term     = c("phi", "theta", "co"),
    estimate = c(x$phi, x$theta, x$co)
  )
}

#' @importFrom fabletools report
#' @rdname NNARMA
#' @export
report.NNARMA <- function(object, ...) {
  cat("  ARMA parameters:\n")
  cat(sprintf("    phi   = %g\n", object$phi))
  cat(sprintf("    theta = %g\n", object$theta))
  cat(sprintf("    co    = %g\n", object$co))
  cat(sprintf("\n  Error variance: %g\n", var(object$v_state)))
  if (!is.null(object$seasons))
    cat(sprintf("  Seasonal period:  %d\n", object$frequency))
  invisible(object)
}

#' Generate sample paths from a NNARMA model
#'
#' @param x A fitted `NNARMA` model object.
#' @inheritParams forecast.NNARMA
#'
#' @return A vector of future paths from a dataset using a fitted model.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, NNARMA(value))
#' generate(fit, new_data = tsibble::new_data(ts, 7))
#' @export
generate.NNARMA <- function(x, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  sim <- unlist(generate(forecast(x, new_data), 1))
  new_data$.sim <- ifelse(sim >= 0, sim, 0)
  new_data
}

nnarma_optimize <- function(y){
  
  # Define the negative log-likelihood function to be optimised
  nnarma_nll <- function(x, y) {
    phi <- x[1]
    theta <- x[2]
    co <- x[3]
    
    process <- armaDynamic(y, phi, theta, co)
    v <- process$v
    sum(v^2)
  }
  
  init_params <- c(1.5, -0.5, 1.0)
  lb <- c(.NNARMA_EPSILON, -Inf, .NNARMA_EPSILON)
  ub <- c(Inf, .NNARMA_EPSILON, Inf)
  opt <- nloptr(
      x0 = init_params,
      eval_f = function(x) nnarma_nll(x, y),
      lb = lb,
      ub = ub,
      eval_g_ineq = function(x) - x[1] - x[2] + .NNARMA_EPSILON,
      opts = list(algorithm = "NLOPT_LN_COBYLA", maxeval = 500)
    )
  opt
}

nnarma_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by NNARMA.")
}