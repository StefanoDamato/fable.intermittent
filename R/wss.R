#' Willemain--Smart--Schwarz Bootstrap Method
#'
#' Markov-chain bootstrap forecasting method for intermittent demand proposed by
#' Willemain, Smart & Schwarz (2004). Occurrence of non-zero demand is modelled
#' with a two-state Markov chain whose transition probabilities are estimated
#' from the observed binary occurrence sequence. Demand sizes are resampled from
#' past non-zero observations with Gaussian jittering to smooth the empirical
#' distribution.
#'
#' @param formula Model specification.
#' @param ... Not used.
#'
#' @references
#' Willemain, T. R., Smart, C. N., & Schwarz, H. F. (2004). A new approach to
#' forecasting intermittent demand for service parts inventories.
#' *International Journal of Forecasting*, 20(3), 375--387.
#' \doi{10.1016/S0169-2070(03)00013-X}.
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
#'   model(WSS(value)) |>
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
#' @importFrom stats rnorm runif
#' @export
WSS <- function(formula, ...) {
  wss_model <- new_model_class(
    "WSS",
    train = train_wss,
    specials = new_specials(
      xreg = no_xreg
    )
  )
  new_model_definition(wss_model, {{ formula }}, ...)
}

train_wss <- function(.data, specials, ...) {

  # Extrapolate values from the tsibble
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by WSS.")
  }
  y <- unclass(.data)[[measured_vars(.data)]]

  # Check for missing values and all-zero series
  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by WSS.")
  }
  if (all(y == 0)) {
    abort("The time series is all zero.")
  }

  # Perform Croston's decomposition
  decomp <- crostons_decomp(y)
  occurrence <- decomp$occurrence
  demand <- decomp$demand

  # Estimate transition probabilities for the Markov chain
  p <- wss_transition_matrix(occurrence)
  last_occurrence <- occurrence[length(occurrence)]

  # Get fitted values and residuals using the Markov chain
  mean_demand <- mean(demand)
  n <- length(y)
  fitted <- rep(NA_real_, n)
  if (n >= 2) {
    prob_occ <- p[cbind(occurrence[1:(n - 1)] + 1, rep(2, n - 1))]
    fitted[2:n] <- prob_occ * mean_demand
  }
  residuals <- y - fitted

  # Save model components in a structured object
  structure(
    list(
      p = p,
      demand = demand,
      last_occurrence = last_occurrence,
      mean_demand = mean_demand,
      fitted = fitted,
      residuals = residuals
    ),
    class = "WSS"
  )
}

#' Forecast a WSS model
#'
#' Produces forecast distributions from a fitted WSS model using simulation.
#'
#' @inheritParams forecast.EMPDISTR
#' @param times The number of sample paths to use in estimating the forecast
#'   distribution.
#'
#' @return A distribution vector of forecasts of class `dist_sample`.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, WSS(value))
#' forecast(fit, h = "7 days")
#'
#' @export
forecast.WSS <- function(object, new_data, specials = NULL, times = 10000, ...) {
  h <- nrow(new_data)
  if (!is_integerish(times) || times <= 0) {
    abort("`times` must be a positive integer.")
  }

  sim <- wss_simulate(object, h, times)
  samples <- as.list(as.data.frame(sim))
  dist_sample(samples)
}

#' Generate sample paths from a WSS model
#'
#' @param x A fitted `WSS` model object.
#' @inheritParams forecast.WSS
#'
#' @return A vector of future paths from a dataset using a fitted model.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, WSS(value))
#' generate(fit, new_data = tsibble::new_data(ts, 7))
#' @export
generate.WSS <- function(x, new_data, specials = NULL, ...) {
  h <- NROW(new_data)
  sim <- wss_simulate(x, h, 1L)
  new_data$.sim <- as.numeric(sim[1, ])
  new_data
}

#' Extract fitted values from a WSS model
#'
#' @inherit fitted.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, WSS(value))
#' fitted(fit)
#' @export
fitted.WSS <- function(object, ...) {
  object$fitted
}

#' Extract residuals from a WSS model
#'
#' @inherit residuals.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, WSS(value))
#' residuals(fit)
#' @export
residuals.WSS <- function(object, ...) {
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
#' fit <- model(ts, WSS(value))
#' model_sum(fit[[1]][[1]])
#' @export
model_sum.WSS <- function(x) {
  "WSS"
}

wss_simulate <- function(object, h, times) {
  forecast_samples <- matrix(0, nrow = times, ncol = h)

  # Set the occurrence state and sample a value
  occ_state <- rep(object$last_occurrence, times)
  for (i in seq_len(h)) {
    p_x1 <- object$p[cbind(occ_state + 1, rep(2, times))]

    # Update occurrence state based on transition probabilities
    occ_state <- as.integer(runif(times) <= p_x1)
    forecast_samples[, i] <- occ_state
  }

  # Replace positive demand via sampling with jittering
  to_sample <- forecast_samples == 1
  if (any(to_sample)) {
    new_demand <- sample(object$demand, sum(to_sample), replace = TRUE)
    jittered <- 1 + floor(new_demand + rnorm(sum(to_sample)) * sqrt(new_demand))
    jittered[jittered <= 0] <- new_demand[jittered <= 0]
    forecast_samples[to_sample] <- jittered
  }

  forecast_samples
}


wss_transition_matrix <- function(occurrence) {
  len <- length(occurrence)
  occ_diff <- 2 * occurrence[2:len] - occurrence[1:(len - 1)]
  p <- matrix(c(
    sum(occ_diff == 0), sum(occ_diff == -1),
    sum(occ_diff == 2), sum(occ_diff == 1)
  ), 2, 2)

  for (i in 1:2) {
    if (sum(p[i, ]) == 0) {
      p[i, ] <- c(length(occurrence) - sum(occurrence), sum(occurrence))
    }
  }

  p / rowSums(p)
}


no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by WSS.")
}
