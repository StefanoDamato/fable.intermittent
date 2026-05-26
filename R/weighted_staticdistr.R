#' Weighted Static Parametric Distribution
#'
#' Parametric count-distribution forecasting method that weights recent
#' observations more heavily when estimating the distribution parameters.
#' Only the last half of the (hot-start-trimmed) series is used.  Within that
#' pool a weight proportional to \eqn{\alpha^{T-t}} is assigned to observation
#' \eqn{y_t}, where \eqn{T} is the pool length and
#' \eqn{\alpha \in (0, 1]} is the decay parameter.
#' Distribution parameters are estimated via weighted method of moments.
#'
#' When \code{alpha} is not supplied it is selected automatically via
#' time-series leave-one-out cross-validation (LOOCV) on the log-score.
#'
#' @section Weighted method of moments:
#' Let \eqn{w_k = \alpha^{T-k}} (unnormalised) for \eqn{k = 1, \ldots, T}.
#' \describe{
#'   \item{Poisson}{\eqn{\hat\lambda = \sum w_k y_k / \sum w_k}.}
#'   \item{Negative binomial}{Weighted mean \eqn{\hat\mu} and biased weighted
#'     variance \eqn{\hat\sigma^2 = \sum w_k y_k^2 / \sum w_k - \hat\mu^2}.
#'     If \eqn{\hat\sigma^2 \le \hat\mu} (underdispersed), falls back to
#'     Poisson.}
#' }
#'
#' @section LOOCV details:
#' For each evaluation point \eqn{t = 2, \ldots, n} in the reference pool,
#' parameters are fitted from \eqn{y_1, \ldots, y_{t-1}} with weight
#' \eqn{\alpha^{t-s}} for \eqn{y_s}.  Running weighted sums are maintained
#' recursively (\eqn{O(n)} per optimizer call).  The log-score at \eqn{y_t}
#' under the refitted distribution is accumulated and maximised over
#' \eqn{\alpha \in (0, 1]} via \code{\link[stats]{optimize}}.
#'
#' @section \code{generate} vs \code{forecast}:
#' \describe{
#'   \item{\code{forecast}}{Returns the distribution fitted from the full
#'     reference pool (identical for every horizon after weight
#'     normalisation).}
#'   \item{\code{generate}}{Recursive: \eqn{y_{T+1}^*} is drawn from the
#'     distribution fitted on \eqn{(y_1,\ldots,y_T)}; at each subsequent
#'     step the drawn value is appended to the pool with the highest weight
#'     and the distribution is refitted before drawing the next value.}
#' }
#'
#' @param formula Model specification.
#' @param distr Distribution family: \code{"pois"} (Poisson) or
#'   \code{"nbinom"} (negative binomial).
#' @param hot_start Logical.  If \code{TRUE}, leading zeros are removed before
#'   fitting.
#' @param alpha Numeric in \eqn{(0, 1]} or \code{NULL} (default).  Decay
#'   parameter.  When \code{NULL}, selected automatically via LOOCV.
#' @param ... Not used.
#'
#' @return A model specification.
#'
#' @importFrom fabletools new_model_class new_specials new_model_definition
#' @importFrom tsibble measured_vars
#' @importFrom rlang abort arg_match
#' @importFrom distributional dist_poisson dist_negative_binomial generate
#' @importFrom stats dpois dnbinom optimize
#' @export
WSTATICDISTR <- function(formula,
                         distr     = c("pois", "nbinom"),
                         hot_start = FALSE,
                         alpha     = NULL,
                         ...) {
  distr <- arg_match(distr)
  if (!is.null(alpha) &&
      (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha > 1)) {
    abort("`alpha` must be a single numeric value in (0, 1] or NULL.")
  }
  wstaticdistr_model <- new_model_class(
    "WSTATICDISTR",
    train    = train_wstaticdistr,
    specials = new_specials(
      xreg = wstaticdistr_no_xreg
    )
  )
  new_model_definition(wstaticdistr_model, {{ formula }},
                       distr = distr, hot_start = hot_start, alpha = alpha, ...)
}


train_wstaticdistr <- function(.data, specials,
                               distr, hot_start = FALSE, alpha = NULL, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by WSTATICDISTR.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by WSTATICDISTR.")
  }

  T_full    <- length(y)
  start     <- if (hot_start) min(which(y != 0)) else 1L
  y_trimmed <- y[start:T_full]

  # Reference pool: last half of the (trimmed) series
  n_trimmed  <- length(y_trimmed)
  half_start <- ceiling(n_trimmed / 2)
  y_emp      <- y_trimmed[half_start:n_trimmed]

  alpha_fit <- if (is.null(alpha)) wstaticdistr_loocv(y_emp, distr) else alpha

  n     <- length(y_emp)
  raw_w <- alpha_fit ^ seq.int(n, 1L, by = -1L)   # y_emp[n] (y_T) gets alpha^1
  fitted_distr <- wstaticdistr_mom_fit(y_emp, raw_w, distr)

  fitted_vals <- rep(mean(fitted_distr), T_full)
  residuals   <- y - fitted_vals

  structure(
    list(
      y_emp        = y_emp,
      alpha        = alpha_fit,
      distr        = distr,
      fitted_distr = fitted_distr,
      T_full       = T_full,
      fitted       = fitted_vals,
      residuals    = residuals
    ),
    class = "WSTATICDISTR"
  )
}


# Weighted method-of-moments fit.
# raw_w can be unnormalised; normalisation is done internally.
wstaticdistr_mom_fit <- function(y, raw_w, distr) {
  Z  <- sum(raw_w)
  mu <- sum(raw_w * y) / Z
  mu <- max(mu, wstaticdistr_epsilon)

  if (distr == "pois") {
    return(dist_poisson(mu))
  }

  # NB: biased weighted variance
  if (length(y) < 2L) {
    return(dist_poisson(mu))
  }
  var_w <- sum(raw_w * y^2) / Z - mu^2

  if (!is.finite(var_w) || var_w <= mu) {
    return(dist_poisson(mu))
  }

  size <- mu^2 / (var_w - mu)
  prob <- mu / var_w
  size <- max(size, wstaticdistr_epsilon)
  prob <- min(max(prob, wstaticdistr_epsilon), 1 - wstaticdistr_epsilon)
  dist_negative_binomial(size, prob)
}


# Time-series LOO cross-validation for the decay parameter alpha.
#
# For t = 2, ..., n the distribution is fitted from y_emp[1:(t-1)] with weight
# alpha^(t-s) for position s.  Running weighted sufficient statistics are
# maintained recursively (O(n) per optimizer call):
#
#   S_t  = sum_{s=1}^{t-1} alpha^(t-s) * y_s          (weighted sum)
#   Q_t  = sum_{s=1}^{t-1} alpha^(t-s) * y_s^2        (weighted sum of squares)
#   Z_t  = sum_{s=1}^{t-1} alpha^(t-s)                 (normalisation)
#
# Recursion (initialised at t=2 with y_emp[1]):
#   S_{t+1} = alpha * (S_t + y_t),  Q_{t+1} = alpha * (Q_t + y_t^2),
#   Z_{t+1} = alpha * (Z_t + 1)
wstaticdistr_loocv <- function(y_emp, distr) {
  n <- length(y_emp)
  if (n < 3L) return(0.9)

  neg_log_score <- function(alpha) {
    # Initialise running stats for evaluating at t = 2 (history = {y_emp[1]})
    Stilde <- alpha * y_emp[1L]
    Qtilde <- alpha * y_emp[1L]^2
    Ztilde <- alpha

    total <- 0
    valid <- 0L

    for (t in 2L:n) {
      mu_loo <- max(Stilde / Ztilde, wstaticdistr_epsilon)

      ld <- if (distr == "pois" || t < 3L) {
        # t < 3: history has 1 point → variance unreliable, use Poisson
        dpois(y_emp[t], mu_loo, log = TRUE)
      } else {
        var_loo <- Qtilde / Ztilde - mu_loo^2
        if (!is.finite(var_loo) || var_loo <= mu_loo) {
          dpois(y_emp[t], mu_loo, log = TRUE)
        } else {
          size_loo <- max(mu_loo^2 / (var_loo - mu_loo), wstaticdistr_epsilon)
          prob_loo <- min(max(mu_loo / var_loo, wstaticdistr_epsilon),
                         1 - wstaticdistr_epsilon)
          dnbinom(y_emp[t], size_loo, prob_loo, log = TRUE)
        }
      }

      if (is.finite(ld)) {
        total <- total + ld
        valid <- valid + 1L
      }

      # Update running stats: add y_emp[t] to the pool for the next step
      Stilde <- alpha * (Stilde + y_emp[t])
      Qtilde <- alpha * (Qtilde + y_emp[t]^2)
      Ztilde <- alpha * (Ztilde + 1)
    }

    if (valid == 0L) return(Inf)
    -total
  }

  opt <- optimize(neg_log_score, interval = c(1e-6, 1 - 1e-6), tol = 1e-4)
  min(max(opt$minimum, 1e-6), 1)
}


#' @export
forecast.WSTATICDISTR <- function(object, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  # The normalised weights are invariant to the alpha^h scaling factor, so
  # the forecast distribution is the same for every horizon.
  rep(object$fitted_distr, h)
}


#' @export
generate.WSTATICDISTR <- function(x, new_data, specials = NULL, ...) {
  h     <- nrow(new_data)
  alpha <- x$alpha
  distr <- x$distr

  # Recursive: the drawn value is appended to the pool so that it becomes
  # the most-recent observation (weight alpha^1) at the next step.
  pool <- x$y_emp
  sim  <- numeric(h)

  for (i in seq_len(h)) {
    m     <- length(pool)
    raw_w <- alpha ^ seq.int(m, 1L, by = -1L)   # pool[m] (newest) gets alpha^1
    d     <- wstaticdistr_mom_fit(pool, raw_w, distr)
    sim[i] <- distributional::generate(d, 1L)[[1L]]
    pool   <- c(pool, sim[i])
  }

  new_data$.sim <- as.numeric(sim)
  new_data
}


#' @export
fitted.WSTATICDISTR <- function(object, ...) {
  object$fitted
}


#' @export
residuals.WSTATICDISTR <- function(object, ...) {
  object$residuals
}


#' @export
model_sum.WSTATICDISTR <- function(x) {
  sprintf("WSTATICDISTR(%s, alpha=%.3f)", x$distr, x$alpha)
}


wstaticdistr_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by WSTATICDISTR.")
}

wstaticdistr_epsilon <- 1e-4
