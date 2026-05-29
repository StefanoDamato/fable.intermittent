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
#' 15(4), 267--280. \doi{10.1007/s11222-005-4070-y}.
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
#' distributional::variance(d)
#' distributional::generate(d, 5)
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

#' @importFrom stats density
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

#' @importFrom distributional generate
#' @exportS3Method distributional::generate
generate.dist_tweedie <- function(x, times, ...) {
  rtweedie(times,
    mean = x[["mu"]],
    dispersion = x[["phi"]],
    power = x[["p"]]
  )
}

#' @exportS3Method distributional::cdf
cdf.dist_tweedie <- function(x, at, lower.tail = TRUE, log.p = FALSE, ...) {
  ptweedie(at,
    mean = x[["mu"]],
    dispersion = x[["phi"]],
    power = x[["p"]],
    lower.tail = lower.tail,
    log.p = log.p
  )
}

#' @exportS3Method distributional::quantile
quantile.dist_tweedie <- function(x, p, lower.tail = TRUE, log.p = FALSE, ...) {
  qtweedie(p,
    mean = x[["mu"]],
    dispersion = x[["phi"]],
    power = x[["p"]],
    lower.tail = lower.tail,
    log.p = log.p
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


#' Tweedie Distribution Functions
#'
#' @description
#' Density, distribution function, quantile function and random generation for
#' the Tweedie distribution with mean equal to `mean`, dispersion equal to
#' `dispersion`, and power equal to `power`.
#'
#' @details
#' If `mean`, `dispersion`, or `power` are not specified they assume the
#' default values of `1`, `1`, and `1.5`, respectively.
#'
#' The Tweedie distribution used here follows the compound Poisson-Gamma
#' parameterisation with power parameter in \eqn{(1, 2)}. It has
#' \eqn{\mathbb{E}[X] = \mu} and
#' \eqn{\mathrm{Var}(X) = \phi\mu^p}, where \eqn{\mu} is `mean`,
#' \eqn{\phi} is `dispersion`, and \eqn{p} is `power`.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken
#'   to be the number required.
#' @param mean vector of means.
#' @param dispersion vector of dispersion parameters.
#' @param power vector of power parameters.
#' @param log,log.p logical; if `TRUE`, probabilities `p` are given as `log(p)`.
#' @param lower.tail logical; if `TRUE` (default), probabilities are
#'   \eqn{P[X \le x]}; otherwise, \eqn{P[X > x]}.
#'
#' @return 
#' `dtweedie` gives the density, `ptweedie` gives the distribution
#'  function, `qtweedie` gives the quantile function, and `rtweedie`
#'  generates random samples.
#' 
#' The length of the result is determined by `n` for `rtweedie`, and is the
#' maximum of the lengths of the numerical arguments for the other functions.
#'
#' The numerical arguments other than `n` are recycled to the length of the
#' result. Only the first elements of the logical arguments are used.
#'
#' @references
#' Dunn, P. K., & Smyth, G. K. (2005). Series evaluation of Tweedie
#' exponential dispersion model densities. *Statistics and Computing*,
#' 15(4), 267--280. \doi{10.1007/s11222-005-4070-y}.
#'
#' @name tweedie
#' @rdname tweedie
#' @aliases dtweedie ptweedie qtweedie rtweedie
#' @export
rtweedie <- function(n, mean = 1, dispersion = 1, power = 1.5) {
  lambda <- (mean^(2 - power)) / (dispersion * (2 - power))
  alpha <- (2 - power) / (power - 1)
  beta <- 1 / (dispersion * (power - 1) * (mean^(power - 1)))

  m <- rpois(n, lambda)
  rgamma(n, m * alpha, beta)
}

#' @rdname tweedie
#' @export
dtweedie <- function(x, mean = 1, dispersion = 1, power = 1.5, log = FALSE) {
  as.vector(tweedieDensity(x, mean, dispersion, power, log))
}

#' @rdname tweedie
#' @export
ptweedie <- function(q, mean = 1, dispersion = 1, power = 1.5, lower.tail = TRUE, log.p = FALSE) {
  cdf <- as.vector(tweedieCDF(q, mean, dispersion, power))
  if (!lower.tail) {
    cdf <- 1 - cdf
  }
  if (log.p) {
    cdf <- log(cdf)
  }
  cdf
}

#' @rdname tweedie
#' @export
qtweedie <- function(p, mean = 1, dispersion = 1, power = 1.5, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  invcdf <- as.vector(tweedieInvCDF(p, mean, dispersion, power))
  invcdf
}
