#' Tweedie Distribution
#'
#' Construct a Tweedie distribution object using the compound Poisson--Gamma
#' parameterisation with power parameter in \eqn{(1, 2)}. The Tweedie family
#' is a subclass of exponential dispersion models that naturally produces exact
#' zeros (via the Poisson count component) mixed with continuous positive
#' values (via the Gamma severity component), making it well suited to
#' intermittent demand data.
#'
#' The density is evaluated using the series expansion of Dunn & Smyth (2005),
#' implemented in C++ for performance.
#'
#' @param mean Mean parameter \eqn{\mu > 0}.
#' @param dispersion Dispersion parameter \eqn{\phi > 0}.
#' @param power Power parameter \eqn{p \in (1, 2)}.
#'
#' @return A `distributional` distribution object of class `dist_tweedie`.
#'
#' @references
#' Dunn, P. K., & Smyth, G. K. (2005). Series evaluation of Tweedie
#' exponential dispersion model densities. *Statistics and Computing*,
#' 15(4), 267--280.
#'
#' @export
#'
#' @importFrom rlang abort
#' @importFrom distributional new_dist covariance
#' @importFrom stats rpois rgamma
#'
#' @examples
#' d <- dist_tweedie(mean = 2, dispersion = 0.8, power = 1.5)
#' mean(d)
#' variance(d)
#' generate(d, 5)
dist_tweedie <- function(mean = 1, dispersion = 1, power = 1.5) {
  mean <- as.double(mean)
  dispersion <- as.double(dispersion)
  power <- as.double(power)

  if (any(mean <= 0, na.rm = TRUE)) {
    abort("The mean parameter of a Tweedie distribution must be strictly positive.")
  }
  if (any(dispersion <= 0, na.rm = TRUE)) {
    abort("The dispersion parameter of a Tweedie distribution must be strictly positive.")
  }
  if (any(power <= 1 | power >= 2, na.rm = TRUE)) {
    abort("The power parameter of a Tweedie distribution must be in (1, 2).")
  }

  new_dist(mu = mean, phi = dispersion, p = power, class = "dist_tweedie")
}

#' @export
format.dist_tweedie <- function(x, digits = 2, ...) {
  sprintf(
    "Tweedie(%s, %s, %s)",
    format(x[["mu"]], digits = digits, ...),
    format(x[["phi"]], digits = digits, ...),
    format(x[["p"]], digits = digits, ...)
  )
}

#' @exportS3Method distributional::density
#' @export
density.dist_tweedie <- function(x, at, ...) {
  dtweedie(at,
    mean = x[["mu"]],
    dispersion = x[["phi"]],
    power = x[["p"]],
    log = FALSE
  )
}

#' @exportS3Method distributional::log_density
#' @export
log_density.dist_tweedie <- function(x, at, ...) {
  dtweedie(at,
    mean = x[["mu"]],
    dispersion = x[["phi"]],
    power = x[["p"]],
    log = TRUE
  )
}

#' @exportS3Method distributional::generate
generate.dist_tweedie <- function(x, times, ...) {
  rtweedie(times,
    mean = x[["mu"]],
    dispersion = x[["phi"]],
    power = x[["p"]]
  )
}

#' @export
mean.dist_tweedie <- function(x, ...) {
  x[["mu"]]
}

#' @export
covariance.dist_tweedie <- function(x, ...) {
  x[["phi"]] * x[["mu"]]^x[["p"]]
}

#' Sample from Tweedie Distribution
#'
#' Uses the compound Poisson--Gamma representation to draw exact samples
#' for power in \eqn{(1, 2)} (Dunn & Smyth, 2005).
#'
#' @param n The number of samples.
#' @param mean Mean parameter of the distribution.
#' @param dispersion Dispersion parameter of the distribution.
#' @param power Power parameter of the distribution.
#'
#' @keywords internal
rtweedie <- function(n, mean = 1, dispersion = 1, power = 1.5) {
  lambda <- (mean^(2 - power)) / (dispersion * (2 - power))
  alpha <- (2 - power) / (power - 1)
  beta <- 1 / (dispersion * (power - 1) * (mean^(power - 1)))

  m <- rpois(n, lambda)
  rgamma(n, m * alpha, beta)
}

#' Compute Tweedie Density
#'
#' Evaluates the Tweedie density via the series expansion of Dunn & Smyth
#' (2005), dispatching to the C++ implementation `tweedieDensity`.
#'
#' @param x The values the density is evaluated at.
#' @param mean Mean parameter of the distribution.
#' @param dispersion Dispersion parameter of the distribution.
#' @param power Power parameter of the distribution.
#' @param log Whether to return the logarithm.
#'
#' @keywords internal
dtweedie <- function(x, mean = 1, dispersion = 1, power = 1.5, log = FALSE) {
  as.vector(tweedieDensity(x, mean, dispersion, power, log))
}
