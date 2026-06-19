#' Markov Chain Model with Random Walk dynamic
#'
#' A Random Walk model on the demand composed by a switching constant and an ARMA model,
#' where the changes are controlled by a Markov Chain on the occurrence process.
#' Parameters are estimated in closed-form on deseasonalized data,
#' and forecasts are returned as Gaussian distributions for each time step.
#'
#' @param formula Model specification.
#' @param ... Not used.
#'
#' @references
#'
#' Sbrana, G. (2025). Markov Walk and Walmart sales prediction.
#' *Journal of the Operational Research Society*, 1--12. \doi{10.1080/01605682.2025.2569661}.
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
#'   model(MARWAL(value)) |>
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
#' @importFrom distributional dist_normal dist_truncated
#' @export
MARWAL <- function(formula, ...) {
  marwal_model <- new_model_class(
    "MARWAL",
    train = train_marwal,
    specials = new_specials(
      xreg = marwal_no_xreg
    )
  )
  new_model_definition(marwal_model, {{ formula }}, ...)
}


train_marwal <-function(.data, specials, ...) {
  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by MARWAL.")
  }
  if (all(y == 0)) {
    abort("The time series is all zero.")
  }

  # Extract the correct frequency and deseasonalize the data
  max_prop_zeros <- 0.95 #TODO: can this be made as an usable parameter?
  period <- get_freq(.data)
  if (period < 1) {
    abort("The seasonal period must be greater than or equal to 1.")
  }
  deseasonalized <- marwal_deseasonalize(y, period = period,
                                         max_prop_zeros = max_prop_zeros)
  y_deseasonalized <- deseasonalized$y_deseasonalized
  seasons <- deseasonalized$seasons

  # Perform Croston's decomposition
  decomp <- crostons_decomp(y_deseasonalized)
  occurrence <- decomp$occurrence
  demand <- decomp$demand

  # Estimate transition probabilities for the Markov chain
  mc <- marwal_transition_matrix(occurrence)
  lambda <- mc$lambda
  xi <- mc$xi
  delta <- mc$delta

  # Estimate the filtering variables
  z <- v <- numeric(length(y))
  mean_demand <- mean(demand)
  k <- (lambda^2 - 1 + sqrt(1 - lambda^2)) / lambda
  for (t in 1:length(y)) {
    v[t] <- y_deseasonalized[t] - occurrence[t] * mean_demand - ifelse(t == 1, 0, z[t - 1])
    z[t] <- lambda * ifelse(t == 1, 0, z[t - 1]) + k * v[t]
  }
  var_v <- var(v)

  # Compute the error variance
  p <- mean(occurrence)
  sigmasq <- ifelse(p < 1, ((-1 + lambda) * (1 + lambda - 2 * p) * p) / (-1 + p), 0)

  # Deseasonalize fitted values and residuals
  fitted <- z * ifelse(is.null(seasons), 1, rep(seasons, length.out = length(y)))
  residuals <- y - fitted

  structure(
    list(
      lambda = lambda,
      xi = xi,
      delta = delta,
      k = k,
      mean_demand = mean_demand,
      var_v = var_v,
      sigmasq = sigmasq,
      frequency = period,
      seasons = if (period > 1) seasons else NULL,
      last_z = z[length(z)],
      last_occurrence = occurrence[length(occurrence)],
      length_y = length(y),
      fitted = fitted,
      residuals = residuals
    ),
    class = "MARWAL"
  )
}


marwal_transition_matrix <- function(occurrence, mean_y) {

  # Compute the transition counts for the Markov chain
  len <- length(occurrence)
  occ_diff <- 2 * occurrence[2:len] - occurrence[1:(len - 1)]
  num <- matrix(c(
    sum(occ_diff == 0), sum(occ_diff == -1),
    sum(occ_diff == 2), sum(occ_diff == 1)
  ), 2, 2)

  # Extract transition probabilities
  lambda <- max(num[2, 2] / sum(num[2, ]), .MARWAL_EPSILON)
  xi <- num[1, 2] / sum(num[1, ])
  p00 <- 1 - xi
  if (sum(num[1, ]) == 0) lambda <- 1
  if (sum(num[2, ]) == 0) lambda <- 0
  delta <- p00 + lambda - 1
  if (lambda == 1) {
    xi <- 0
    delta <- 1
  }

  list(lambda = lambda, xi = xi, delta = delta)
}


marwal_deseasonalize <- function(y, period, max_prop_zeros) {

  # Ignore the seasonality if there are too many zeros
  if ((period <= 1) | ((length(y[y == 0]) / length(y)) >= max_prop_zeros)) {
    return(list(y_deseasonalized = y, seasons = NULL))
  }

  # Compute residuals of the moving average
  moving_avg <- rep(NA, length(y))
  for (i in 1:(length(y) - period + 1)) {
    moving_avg[i + ((period + 1) / 2) - 1] <- mean(y[i:(i + period - 1)])
  }
  resid <- ifelse(moving_avg > 0, y / moving_avg, 0)

  # Compute the seasonal factors and deseasonalise the data
  seasons <- numeric(period)
  for (s in 1:period) {
    seasons[s] <- mean(resid[seq(s, length(y) - period + s, by = period)], na.rm = TRUE)
  }
  y_deseasonalized <- y / rep(seasons, length.out = length(y))
  y_deseasonalized[!is.finite(y_deseasonalized)] <- 0

  list(y_deseasonalized = y_deseasonalized, seasons = seasons)
}

#' Forecast a MARWAL model
#'
#' Produces forecast distributions from a fitted MARWAL model.
#'
#' @inheritParams forecast.EMPDISTR
#'
#' @return A distribution vector of class `dist_normal`.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, MARWAL(value))
#' forecast(fit, h = "7 days")
#'
#' @export
forecast.MARWAL <- function(object, new_data, specials = NULL, ...) {
  h <- nrow(new_data)

  # Compute the mean forecast
  mc_state <- object$last_occurrence

  # Compute the mean forecasts
  mc_fc <- (object$delta^(1:h) * object$last_occurrence +
              object$xi * cumsum(object$delta^(0:(h - 1))))
  mean_fc <- mc_fc * object$mean_demand + object$lambda^(0:(h - 1)) * object$last_z


  # Adjust for seasonality if necessary
  if (!is.null(object$seasons)) {
    s <- 1 + object$length_y %% object$frequency
    seasons <- rep(object$seasons,
                   1 + ceiling((h - object$frequency + s - 1) / object$frequency))
    mean_fc <- mean_fc * seasons[1:h + s -1]
  }

  # Calculate forecast variance
  delta_sum <- c(1, cumsum(object$delta^(2 * (0:(h - 2)))))
  lambda_sum <- c(1, cumsum(object$lambda^(2 * (0:(h - 2)))))
  var_fc <- object$sigmasq * delta_sum * object$mean_demand^2 +
    object$var_v * (1 + object$k^2 * lambda_sum)

  # Return the Gaussian forecast distribution
  dist_truncated(dist_normal(mean_fc, sqrt(var_fc)), lower = 0)
}

#' Extract fitted values from a MARWAL model
#'
#' @inherit fitted.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, MARWAL(value))
#' fitted(fit)
#'
#' @export
fitted.MARWAL <- function(object, ...) {
  object$fitted
}

#' Extract residuals from a MARWAL model
#'
#' @inherit residuals.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, MARWAL(value))
#' residuals(fit)
#' @export
residuals.MARWAL <- function(object, ...) {
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
#' fit <- model(ts, MARWAL(value))
#' model_sum(fit[[1]][[1]])
#' @export
model_sum.MARWAL <- function(x) {
  "MARWAL"
}

#' Generate sample paths from a MARWAL model
#'
#' @param x A fitted `MARWAL` model object.
#' @inheritParams forecast.MARWAL
#'
#' @return A vector of future paths from a dataset using a fitted model.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, MARWAL(value))
#' generate(fit, new_data = tsibble::new_data(ts, 7))
#' @export
generate.MARWAL <- function(x, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  sim <- unlist(generate(forecast(x, new_data), 1))
  new_data$.sim <- ifelse(sim >= 0, sim, 0)
  new_data
}

marwal_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by MARWAL.")
}