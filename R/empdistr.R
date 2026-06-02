#' Empirical Distribution Resampling
#'
#' Naive non-parametric baseline for intermittent demand forecasting. The
#' predictive distribution at every horizon is simply the empirical distribution
#' of the observed values: forecasts are produced by resampling with
#' replacement from the historical series. Point forecasts are the sample mean.
#'
#' @param formula Model specification.
#' @param hot_start Logical. If `TRUE`, leading zeros are removed from the
#'   time series before fitting.
#' @param ... Not used.
#'
#' @references
#' Hasni, M., Aguir, M. S., Babai, M. Z., & Jemai, Z. (2019). Spare parts
#' demand forecasting: a review on bootstrapping methods. *International
#' Journal of Production Research*, 57(15--16), 4791--4804.
#' \doi{10.1080/00207543.2018.1424375}.
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
#'   model(EMPDISTR(value)) |>
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
#' @export
EMPDISTR <- function(formula, hot_start = FALSE, ...) {
  empdistr_model <- new_model_class(
    "EMPDISTR",
    train = train_empdistr,
    specials = new_specials(
      xreg = empdistr_no_xreg
    )
  )
  new_model_definition(empdistr_model, {{ formula }}, hot_start = hot_start, ...)
}

train_empdistr <- function(.data, specials, hot_start = FALSE, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by empdistr.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by empdistr.")
  }

  # Remove leading zeros for hot_start
  start <- ifelse(hot_start, min(which(y != 0)), 1)
  y_emp <- y[start:length(y)]

  # Fit the model by simply repeating the mean
  fitted <- rep(mean(y_emp), length(y))
  residuals <- y - fitted

  structure(
    list(
      y_emp = y_emp,
      fitted = fitted,
      residuals = residuals
    ),
    class = "EMPDISTR"
  )
}

#' Forecast an EMPDISTR model
#'
#' Produces forecast distributions by repeating the empirical distribution
#' estimated from the training data at each forecast horizon.
#'
#' @inheritParams generics::forecast
#' @param new_data A tsibble containing future index values to forecast.
#' @param specials Passed by [fabletools::forecast.mdl_df()].
#'
#' @return A distribution vector of class `dist_sample`.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, EMPDISTR(value))
#' forecast(fit, h = "7 days")
#' @export
forecast.EMPDISTR <- function(object, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  samples <- rep(list(object$y_emp), h)
  dist_sample(samples)
}

#' Generate sample paths from an EMPDISTR model
#'
#' @param x A fitted `EMPDISTR` model object.
#' Simulates future observations by resampling with replacement from the
#' empirical support learned during training.
#'
#' @inheritParams forecast.EMPDISTR
#' @param x A fitted `EMPDISTR` model object.
#'
#' @return A `new_data` tibble with a `.sim` column of simulated values.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, EMPDISTR(value))
#' generate(fit, h = 7)
#' @export
generate.EMPDISTR <- function(x, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  sim <- sample(x$y_emp, size = h, replace = TRUE)
  new_data$.sim <- as.numeric(sim)
  new_data
}

#' Extract fitted values from an EMPDISTR model
#'
#' @param object A model for which fitted values are required.
#' @param ... Not used.
#'
#' @return A numeric vector of fitted values.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, EMPDISTR(value))
#' fitted(fit)
#' @export
fitted.EMPDISTR <- function(object, ...) {
  object$fitted
}

#' Extract residuals from an EMPDISTR model
#'
#' @param object A model for which residuals are required.
#' @param ... Not used.
#'
#' @return A numeric vector of residuals.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, EMPDISTR(value))
#' residuals(fit)
#' @export
residuals.EMPDISTR <- function(object, ...) {
  object$residuals
}

#' Return model name
#' @param x A fitted `fable` model object.
#' @return The model name as a string
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, EMPDISTR(value))
#' model_sum(fit[[1]][[1]])
#' 
#' @importFrom fabletools model_sum
#' @export
model_sum.EMPDISTR <- function(x) {
  "EMPDISTR"
}

empdistr_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by EMPDISTR.")
}
