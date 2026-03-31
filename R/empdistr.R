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
#'
#' @return A model specification.
#'
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

#' @export
forecast.EMPDISTR <- function(object, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  samples <- rep(list(object$y_emp), h)
  dist_sample(samples)
}

#' @export
generate.EMPDISTR <- function(x, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  sim <- sample(x$y_emp, size = h, replace = TRUE)
  new_data$.sim <- as.numeric(sim)
  new_data
}

#' @export
fitted.EMPDISTR <- function(object, ...) {
  object$fitted
}

#' @export
residuals.EMPDISTR <- function(object, ...) {
  object$residuals
}

#' @export
model_sum.EMPDISTR <- function(x) {
  "EMPDISTR"
}

empdistr_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by EMPDISTR.")
}
