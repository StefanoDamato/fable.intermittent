################################################################################
# GLOBAL PARAMETERS (EPSILON) TO AVOID NUMERICAL ISSUES IN COMPUTATIONS
#' @importFrom distributional dist_inflated dist_transformed new_dist covariance
#' @importFrom fabletools get_frequencies
#' @importFrom rlang abort
#' @importFrom stats dnorm pnorm qnorm rnorm
NULL

.BETANBB_EPSILON     <- 1e-4
.GAMPOISB_EPSILON    <- 1e-4
.HSPES_EPSILON       <- 1e-4
.MARWAL_EPSILON      <- 1e-4
.NNARMA_EPSILON      <- 1e-4
.NEGBINES_EPSILON    <- 1e-4
.STATICDISTR_EPSILON <- 1e-4
.TWEES_EPSILON       <- 1e-4

crostons_decomp <- function(y) {
  occurrence <- ifelse(y > 0, 1L, 0L)
  d_times <- which(y > 0)
  demand <- y[d_times]
  intervals <- diff(c(0, d_times))

  list(
    occurrence = occurrence,
    demand = demand,
    intervals = intervals
  )
}

get_freq <- function(.data, period = NULL, model_name = "Model") {
  period <- get_frequencies(period, .data)
  period <- round(as.numeric(period[[1]]))

  period <- as.integer(period)
  if (period < 1) {
    abort("The seasonal period must be greater than or equal to 1.")
  }

  period
}

# Classical multiplicative seasonal adjustment via centered moving average,
# shared by NNARMA (normalize = TRUE, as in its reference implementation) and
# MARWAL (normalize = FALSE, matching the Markov Walk reference)
deseasonalize <- function(y, period, max_prop_zeros, normalize = FALSE) {

  # Ignore the seasonality if there are too many zeros
  if ((period <= 1) | ((length(y[y == 0]) / length(y)) >= max_prop_zeros)) {
    return(list(y_deseasonalized = y, seasons = NULL))
  }

  # Compute residuals of the moving average
  moving_avg <- rep(NA, length(y))
  for (i in 1:(length(y) - period + 1)) {
    moving_avg[i + ((period + 1) / 2) - 1] <- mean(y[i:(i + period - 1)])
  }
  # All-zero windows carry no seasonal information: drop them (NA) rather
  # than counting them as ratio 0, as in the reference implementations
  resid <- ifelse(moving_avg > 0, y / moving_avg, NA)

  # Compute the seasonal factors and deseasonalise the data
  seasons <- numeric(period)
  for (s in 1:period) {
    seasons[s] <- mean(resid[seq(s, length(y) - period + s, by = period)], na.rm = TRUE)
  }
  if (normalize) {
    seasons <- seasons * period / sum(seasons)
  }

  # Only deseasonalize when all factors are strictly positive and finite,
  # so training and forecasting always operate on the same scale
  if (!all(is.finite(seasons)) || min(seasons) <= 0) {
    return(list(y_deseasonalized = y, seasons = NULL))
  }

  list(
    y_deseasonalized = y / rep(seasons, length.out = length(y)),
    seasons = seasons
  )
}

make_hurdle_shifted_distr <- function(distr, pzero){
  distr <- dist_transformed(distr, function(x) x + 1, function(x) x - 1)
  dist_inflated(distr, pzero, 0)
}

# Former generic construction of the censored distribution, replaced by
# dist_normal_nonneg below. Unlike the class, its mean() was the coherent
# expectation E[max(X, 0)] rather than the clamped max(mu, 0).
# negative_mass_to_zero <- function(distr) {
#   do.call(c, lapply(as.list(distr), \(d) {
#     p0 <- cdf(d, 0)[[1]]
#     dist_truncated(d, lower = 0) |> dist_inflated(prob = p0)
#   }))
# }

# Gaussian forecast distribution with its negative part collapsed to zero,
# used by NNARMA and MARWAL. Quantiles, cdf, and samples are those of the
# rectified Gaussian max(X, 0); as a deliberate reporting choice, mean() is
# the clamped max(mu, 0) (the distribution's median) rather than the coherent
# expectation E[max(X, 0)], so point forecasts match the clamped Gaussian mean
# of the reference literature. variance() is the Gaussian sigma^2. Hence
# mean() does not match the average of generate() samples when P(X < 0) > 0.
dist_normal_nonneg <- function(mu = 0, sigma = 1) {
  mu <- as.double(mu)
  sigma <- as.double(sigma)

  if (any(sigma <= 0, na.rm = TRUE)) {
    abort("The sigma parameter of a non-negative Gaussian distribution must be strictly positive.")
  }

  new_dist(mu = mu, sigma = sigma, class = "dist_normal_nonneg")
}

#' @noRd
#' @export
format.dist_normal_nonneg <- function(x, digits = 2, ...) {
  sprintf(
    "N+(%s, %s)",
    format(x[["mu"]], digits = digits, ...),
    format(x[["sigma"]]^2, digits = digits, ...)
  )
}

#' @importFrom stats density
#' @exportS3Method distributional::density
#' @export
#' @noRd
density.dist_normal_nonneg <- function(x, at, ...) {
  # At zero the distribution has an atom; report its probability mass there
  ifelse(
    at < 0,
    0,
    ifelse(
      at == 0,
      pnorm(0, x[["mu"]], x[["sigma"]]),
      dnorm(at, x[["mu"]], x[["sigma"]])
    )
  )
}

#' @importFrom distributional generate
#' @exportS3Method distributional::generate
#' @noRd
generate.dist_normal_nonneg <- function(x, times, ...) {
  pmax(rnorm(times, x[["mu"]], x[["sigma"]]), 0)
}

#' @exportS3Method distributional::cdf
#' @noRd
cdf.dist_normal_nonneg <- function(x, q, lower.tail = TRUE, log.p = FALSE, ...) {
  cdf <- ifelse(q < 0, 0, pnorm(q, x[["mu"]], x[["sigma"]]))
  if (!lower.tail) {
    cdf <- 1 - cdf
  }
  if (log.p) {
    cdf <- log(cdf)
  }
  cdf
}

#' @exportS3Method distributional::quantile
#' @noRd
quantile.dist_normal_nonneg <- function(x, p, lower.tail = TRUE, log.p = FALSE, ...) {
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  pmax(qnorm(p, x[["mu"]], x[["sigma"]]), 0)
}

#' @export
#' @noRd
mean.dist_normal_nonneg <- function(x, ...) {
  max(x[["mu"]], 0)
}

#' @export
#' @noRd
covariance.dist_normal_nonneg <- function(x, ...) {
  x[["sigma"]]^2
}

fit_nbinom <- function(y) {
  if (length(y) == 0 || all(y == 0)) {
    return(c(size = 100, prob = 1 - .STATICDISTR_EPSILON))
  }

  fit <- tryCatch(
    nloptr(
      x0 = c(max(mean(y), .STATICDISTR_EPSILON), 0.5),
      eval_f = function(x) -mean(dnbinom(y, x[1], x[2], log = TRUE)),
      lb = c(.STATICDISTR_EPSILON, .STATICDISTR_EPSILON),
      ub = c(Inf, 1 - .STATICDISTR_EPSILON),
      opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || is.null(fit$solution)) {
    mu <- mean(y)
    sigmasq <- var(y)
    if (!is.na(sigmasq) && sigmasq > mu + .STATICDISTR_EPSILON) {
      size <- (mu^2) / (sigmasq - mu)
    } else {
      size <- 100
    }
    prob <- min(size / (size + mu), 1 - .STATICDISTR_EPSILON)
    return(c(size = size, prob = prob))
  }

  c(size = fit$solution[1], prob = fit$solution[2])
}
