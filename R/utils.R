################################################################################
# GLOBAL PARAMETERS (EPSILON) TO AVOID NUMERICAL ISSUES IN COMPUTATIONS
#' @importFrom distributional dist_inflated dist_transformed new_dist covariance
#' @importFrom fabletools get_frequencies
#' @importFrom nloptr nloptr
#' @importFrom rlang abort
#' @importFrom stats dnbinom dnorm pnorm qnorm rnorm var
#' @importFrom tweedieDistr dtweedie ptweedie qtweedie rtweedie
NULL

.BETANBB_EPSILON     <- 1e-4
.GAMPOISB_EPSILON    <- 1e-4
.HSPES_EPSILON       <- 1e-4
.MARWAL_EPSILON      <- 1e-4
.NNARMA_EPSILON      <- 1e-4
.NEGBINES_EPSILON    <- 1e-4
.PARAMSD_EPSILON <- 1e-4
.TWEES_EPSILON       <- 1e-4


.TWEEDIE_POWER_MIN <- 1.2
.TWEEDIE_POWER_MAX <- 1.8

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

#' @importFrom fabletools get_frequencies
#' @importFrom rlang abort
get_freq <- function(.data, period = NULL, model_name = "Model") {
  period <- get_frequencies(period, .data)
  period <- round(as.numeric(period[[1]]))

  period <- as.integer(period)
  if (period < 1) {
    abort("The seasonal period must be greater than or equal to 1.")
  }

  period
}

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
  resid <- ifelse(moving_avg > 0, y / moving_avg, NA)

  # Compute the seasonal factors and deseasonalise the data
  seasons <- numeric(period)
  for (s in 1:period) {
    seasons[s] <- mean(resid[seq(s, length(y) - period + s, by = period)], na.rm = TRUE)
  }
  if (normalize) {
    seasons <- seasons * period / sum(seasons)
  }

  # Only deseasonalize when all factors are strictly positive and finite
  if (!all(is.finite(seasons)) || min(seasons) <= 0) {
    return(list(y_deseasonalized = y, seasons = NULL))
  }

  list(
    y_deseasonalized = y / rep(seasons, length.out = length(y)),
    seasons = seasons
  )
}

#' @importFrom distributional dist_transformed dist_inflated
make_hurdle_shifted_distr <- function(distr, pzero){
  distr <- dist_transformed(distr, function(x) x + 1, function(x) x - 1)
  dist_inflated(distr, pzero, 0)
}


#' @importFrom rlang abort
#' @importFrom distributional new_dist
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

#' @importFrom stats density dnorm pnorm
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
#' @importFrom stats rnorm
#' @exportS3Method distributional::generate
#' @noRd
generate.dist_normal_nonneg <- function(x, times, ...) {
  pmax(rnorm(times, x[["mu"]], x[["sigma"]]), 0)
}

#' @importFrom stats pnorm
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

#' @importFrom stats qnorm
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


#' @importFrom rlang abort
#' @importFrom distributional new_dist
dist_tweedie_discrete <- function(mean = 1, dispersion = 1, power = 1.5) {
  mean <- as.double(mean)
  dispersion <- as.double(dispersion)
  power <- as.double(power)

  if (any(mean <= 0, na.rm = TRUE)) {
    abort("The mean parameter of a discretised Tweedie distribution must be strictly positive.")
  }
  if (any(dispersion <= 0, na.rm = TRUE)) {
    abort("The dispersion parameter of a discretised Tweedie distribution must be strictly positive.")
  }
  if (any(power <= 1 | power >= 2, na.rm = TRUE)) {
    abort("The power parameter of a discretised Tweedie distribution must be in (1, 2).")
  }

  new_dist(mu = mean, phi = dispersion, p = power, class = "dist_tweedie_discrete")
}


#' @importFrom tweedieDistr ptweedie
tweedie_discrete_pmf <- function(k, mu, phi, power) {
  out <- numeric(length(k))
  ok <- is.finite(k) & k >= 0 & k == round(k)
  if (any(ok)) {
    kk <- k[ok]
    upper <- ptweedie(kk + 0.5, mean = mu, dispersion = phi, power = power)
    lower <- ptweedie(pmax(kk - 0.5, 0), mean = mu, dispersion = phi, power = power)
    lower[kk == 0] <- 0
    out[ok] <- pmax(upper - lower, 0)
  }
  out
}


tweedie_discrete_support <- function(mu, phi, power, n_sd = 15, max_k = 1e5) {
  sd <- sqrt(phi * mu^power)
  0:min(max_k, max(10, ceiling(mu + n_sd * sd)))
}

#' @noRd
#' @export
format.dist_tweedie_discrete <- function(x, digits = 2, ...) {
  sprintf(
    "TweedieD(%s, %s, %s)",
    format(x[["mu"]], digits = digits, ...),
    format(x[["phi"]], digits = digits, ...),
    format(x[["p"]], digits = digits, ...)
  )
}

#' @importFrom stats density
#' @exportS3Method distributional::density
#' @export
#' @noRd
density.dist_tweedie_discrete <- function(x, at, ...) {
  tweedie_discrete_pmf(at, x[["mu"]], x[["phi"]], x[["p"]])
}

#' @importFrom distributional generate
#' @importFrom tweedieDistr rtweedie
#' @exportS3Method distributional::generate
#' @noRd
generate.dist_tweedie_discrete <- function(x, times, ...) {
  round(rtweedie(times, mean = x[["mu"]], dispersion = x[["phi"]], power = x[["p"]]))
}

#' @importFrom tweedieDistr ptweedie
#' @exportS3Method distributional::cdf
#' @noRd
cdf.dist_tweedie_discrete <- function(x, q, lower.tail = TRUE, log.p = FALSE, ...) {
  cdf <- ptweedie(pmax(floor(q) + 0.5, 0), mean = x[["mu"]],
                  dispersion = x[["phi"]], power = x[["p"]])
  cdf[q < 0] <- 0
  if (!lower.tail) {
    cdf <- 1 - cdf
  }
  if (log.p) {
    cdf <- log(cdf)
  }
  cdf
}

#' @importFrom tweedieDistr ptweedie qtweedie
#' @exportS3Method distributional::quantile
#' @noRd
quantile.dist_tweedie_discrete <- function(x, p, lower.tail = TRUE, log.p = FALSE, ...) {
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  mu <- x[["mu"]]
  phi <- x[["phi"]]
  power <- x[["p"]]

  step_cdf <- function(k) ptweedie(k + 0.5, mean = mu, dispersion = phi, power = power)

  vapply(p, function(pi) {
    if (is.na(pi)) return(NA_real_)
    if (pi <= 0) return(0)
    if (pi >= 1) return(Inf)
    # Start from the continuous quantile, then step onto the integer lattice
    k <- max(0, round(qtweedie(pi, mean = mu, dispersion = phi, power = power)))
    while (k > 0 && step_cdf(k - 1) >= pi) k <- k - 1
    while (step_cdf(k) < pi) k <- k + 1
    as.double(k)
  }, numeric(1))
}

#' @export
#' @noRd
mean.dist_tweedie_discrete <- function(x, ...) {
  k <- tweedie_discrete_support(x[["mu"]], x[["phi"]], x[["p"]])
  sum(k * tweedie_discrete_pmf(k, x[["mu"]], x[["phi"]], x[["p"]]))
}

#' @export
#' @noRd
covariance.dist_tweedie_discrete <- function(x, ...) {
  k <- tweedie_discrete_support(x[["mu"]], x[["phi"]], x[["p"]])
  pk <- tweedie_discrete_pmf(k, x[["mu"]], x[["phi"]], x[["p"]])
  m <- sum(k * pk)
  sum(k^2 * pk) - m^2
}

#' @importFrom nloptr nloptr
#' @importFrom stats dnbinom var
fit_nbinom <- function(y) {
  if (length(y) == 0 || all(y == 0)) {
    return(c(size = 100, prob = 1 - .PARAMSD_EPSILON))
  }

  fit <- tryCatch(
    nloptr(
      x0 = c(max(mean(y), .PARAMSD_EPSILON), 0.5),
      eval_f = function(x) -mean(dnbinom(y, x[1], x[2], log = TRUE)),
      lb = c(.PARAMSD_EPSILON, .PARAMSD_EPSILON),
      ub = c(Inf, 1 - .PARAMSD_EPSILON),
      opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || is.null(fit$solution)) {
    mu <- mean(y)
    sigmasq <- var(y)
    if (!is.na(sigmasq) && sigmasq > mu + .PARAMSD_EPSILON) {
      size <- (mu^2) / (sigmasq - mu)
    } else {
      size <- 100
    }
    prob <- min(size / (size + mu), 1 - .PARAMSD_EPSILON)
    return(c(size = size, prob = prob))
  }

  c(size = fit$solution[1], prob = fit$solution[2])
}


#' @importFrom rlang abort
#' @importFrom nloptr nloptr
#' @importFrom stats var
#' @importFrom tweedieDistr dtweedie
fit_tweedie <- function(y, discrete = FALSE) {
  if (length(y) == 0 || all(y == 0)) {
    return(c(mean = .PARAMSD_EPSILON, dispersion = 1, power = 1.5))
  }

  mu <- max(mean(y), .PARAMSD_EPSILON)
  phi_start <- max(var(y) / mu, .PARAMSD_EPSILON)

  if (discrete) {
    if (any(y < 0) || any(y != round(y))) {
      abort(paste0(
        "The discretised Tweedie requires non-negative integer observations. ",
        "For a non-count series use `tweedie_discrete = FALSE` together with ",
        "`distr = \"mixture\"` or `distr = \"tweedie\"`."
      ))
    }
    counts <- table(y)
    vals <- as.numeric(names(counts))
    weights <- as.numeric(counts)
    eval_f <- function(x) {
      pmf <- tweedie_discrete_pmf(vals, x[1], x[2], x[3])
      -sum(weights * log(pmax(pmf, .Machine$double.xmin))) / length(y)
    }
    x0 <- c(mu, phi_start, 1.5)
    lb <- c(.PARAMSD_EPSILON, .PARAMSD_EPSILON,
            .TWEEDIE_POWER_MIN + .PARAMSD_EPSILON)
    ub <- c(max(y) * 10 + 1, Inf, .TWEEDIE_POWER_MAX - .PARAMSD_EPSILON)
  } else {
    eval_f <- function(x) {
      -mean(dtweedie(y, mean = mu, dispersion = x[1], power = x[2], log = TRUE))
    }
    x0 <- c(phi_start, 1.5)
    lb <- c(.PARAMSD_EPSILON, .TWEEDIE_POWER_MIN + .PARAMSD_EPSILON)
    ub <- c(Inf, .TWEEDIE_POWER_MAX - .PARAMSD_EPSILON)
  }

  fit <- tryCatch(
    nloptr(
      x0 = x0, eval_f = eval_f, lb = lb, ub = ub,
      opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || is.null(fit$solution)) {
    return(c(mean = mu, dispersion = phi_start, power = 1.5))
  }

  if (discrete) {
    c(mean = fit$solution[1], dispersion = fit$solution[2], power = fit$solution[3])
  } else {
    c(mean = mu, dispersion = fit$solution[1], power = fit$solution[2])
  }
}
