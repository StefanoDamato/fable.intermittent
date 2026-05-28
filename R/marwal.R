#' Markov Chain Model with Random Walk dynamic
#'
#' A Random Walk model on the demand composed by a swithcing constant and an ARMA model,
#' where the changes are controlled by a Markov Chain on the occurrence process.
#' Paramter are estimated in closed-form on deseasonalized data,
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
#' ts |>
#'   model(MARWAL(value)) |>
#'   forecast(h = "7 days") |>
#'   ggtime::autoplot(ts)
#'
#' @importFrom fabletools new_model_class new_specials new_model_definition
#' @importFrom tsibble measured_vars
#' @importFrom rlang abort is_integerish
#' @importFrom distributional dist_normal
#' @export
MARWAL <- function(formula, ...) {
  marwal_model <- new_model_class(
    "MARWAL",
    train = train_marwal,
    specials = new_specials(
      xreg = no_xreg
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

  # S: OLDER VERSION, SLOWER
  #mean_fc <- numeric(h)
  #for (i in seq_len(h)) {
    #mc_state <- object$xi + object$delta * mc_state
    #mean_fc[i] <- object$mean_demand * mc_state + object$lambda^(i - 1) * object$last_z
  #}

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
  dist_normal(mean_fc, sqrt(var_fc))
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



# # ===========================================================================
# Bernoulli <- function(y, steps) {
#   alpha <- -1
#   co <- mean(y[y > 0])

#   bern <- y
#   bern[bern > 0] <- 1

#   su <- bern[-1] + bern[-length(y)]
#   di <- bern[-1] - bern[-length(y)]

#   n00 <- length(su[su == 0])
#   n11 <- length(su[su == 2])
#   n01 <- length(di[di == -1])
#   n10 <- length(di[di == 1])

#   p00 <- n00 / (n00 + n10)
#   xi  <- n10 / (n10 + n00)
#   p   <- mean(bern)
#   lambda <- n11 / (n11 + n01)

#   if (n00 == 0 & n10 == 0) lambda <- 1
#   if (n11 == 0 & n01 == 0) lambda <- 0
#   if (lambda <= 0) lambda <- .0001

#   delta <- p00 + lambda - 1
#   if (lambda == 1) {
#     p <- 1
#     xi <- 0
#     delta <- 1
#   }

#   # =================================
#   # STE: to be added in the forecast method
#   MC <- c()
#   MC[1] <- xi + delta * tail(bern, 1)
#   for (s in 2:steps) {
#     MC[s] <- xi + delta * MC[s - 1]
#   }
#   # =================================

#   m <- v <- c()
#   m[1] <- 0
#   k <- (lambda^2 - 1 + sqrt(1 - lambda^2)) / lambda

#   for (t in 1:length(y)) {
#     v[t] <- y[t] - bern[t] * co - m[t]
#     m[t + 1] <- lambda * m[t] + k * v[t]
#   }

#   # =================================
#   # STE: to be added in the forecast method
#   fo <- c()
#   for (s in 1:steps) {
#     fo[s] <- co * MC[s] + lambda^(s - 1) * m[length(m)]
#   }
#   # =================================

#   list(bern, p, lambda, MC, delta, xi, fo, v, k)
# }

# MW <- function(y, steps) {
#   prob1 <- .5
#   prob2 <- .67
#   prob3 <- .95
#   prob4 <- .99

#   mysum <- function(x, steps) {
#     d <- 0
#     for (f in 0:(steps - 2)) d <- d + x^(f * 2)
#     d
#   }

#   # Classical multiplicative seasonal adjustment
#   if ((length(y[y == 0]) / length(y)) < .95) {
#     s <- 7
#     h <- steps
#     cma <- matrix(NA, length(y), 1)

#     for (g in 1:(length(y) - s + 1)) {
#       cma[g + ((s + 1) / 2) - 1] <- mean(y[g:(g + s - 1)])
#     }

#     residuals <- y / cma

#     sfactors <- c()
#     for (seas in 1:s) {
#       sfactors[seas] <- mean(na.omit(residuals[seq(seas, length(y) - s + seas, by = s)]))
#     }

#     sfactout <- rep(sfactors, length(y) + h)[(length(y) + 1):(length(y) + h)]
#     y <- y / rep(sfactors, ceiling(length(y) / s))[1:length(y)]
#     y[is.na(y)] <- 0
#     y[y == Inf] <- 0

#   } else {
#     sfactout <- rep(1, steps)
#   }

#   h <- c()
#   le <- 14

#   for (i in 0:1000) {
#     Y <- tail(y, (length(y) - i * le))
#     if (length(Y) < 100) break
#     ins <- head(Y, length(Y) - steps)

#     if ((length(ins[ins == 0]) / length(ins)) < .99 & length(ins) > 100) {
#       h[i + 1] <- mean((tail(Y, steps) - Bernoulli(ins, steps)[[7]])^2)
#     }
#   }

#   if (length(h) != 0) {
#     y <- tail(y, (length(y) - which.min(h[!is.na(h)]) * le) + steps)
#   }

#   co <- mean(y[y > 0])
#   ma <- Bernoulli(y, steps)

#   bern   <- ma[[1]]
#   p      <- ma[[2]]
#   lambda <- ma[[3]]
#   MC     <- ma[[4]]
#   delta  <- ma[[5]]
#   fo     <- ma[[7]]
#   v      <- ma[[8]]
#   k      <- ma[[9]]

#   fo <- fo * sfactout

#   if (p < 1) {
#     vari <- ((-1 + lambda) * (1 + lambda - 2 * p) * p) / (-1 + p)
#   } else {
#     vari <- 0
#   }

#   Interv <- c()
#   Interv[1] <- vari * co^2 + var(v)
#   for (j in 2:steps) {
#     Interv[j] <- vari * mysum(delta, j) * co^2 +
#                  (var(v) * (1 + k^2 * (mysum(lambda, j))))
#   }

#   lower0 <- fo
#   lower50 <- fo - qnorm((1 + prob1) / 2) * sqrt(Interv)
#   lower67 <- fo - qnorm((1 + prob2) / 2) * sqrt(Interv)
#   lower95 <- fo - qnorm((1 + prob3) / 2) * sqrt(Interv)
#   lower99 <- fo - qnorm((1 + prob4) / 2) * sqrt(Interv)

#   upper0 <- fo
#   upper50 <- fo + qnorm((1 + prob1) / 2) * sqrt(Interv)
#   upper67 <- fo + qnorm((1 + prob2) / 2) * sqrt(Interv)
#   upper95 <- fo + qnorm((1 + prob3) / 2) * sqrt(Interv)
#   upper99 <- fo + qnorm((1 + prob4) / 2) * sqrt(Interv)

#   list(
#     mean = fo,
#     lower = cbind(lower0, lower50, lower67, lower95, lower99),
#     upper = cbind(upper0, upper50, upper67, upper95, upper99)
#   )
# }
