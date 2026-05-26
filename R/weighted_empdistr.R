#' Weighted Empirical Distribution
#'
#' Non-parametric forecasting method for intermittent demand that weights
#' recent observations more heavily. Only the last half of the (hot-start-
#' trimmed) series is used as the reference pool. Within that pool, the weight
#' assigned to observation \eqn{y_t} is proportional to
#' \eqn{\alpha^{T-t}}, where \eqn{T} is the pool length and
#' \eqn{\alpha \in (0, 1]} is the decay parameter.  Setting \eqn{\alpha = 1}
#' recovers the plain empirical distribution.
#'
#' When \code{alpha} is not supplied it is selected automatically via
#' time-series leave-one-out cross-validation (LOOCV) on the log-score.
#'
#' @section LOOCV details:
#' For each evaluation point \eqn{t = 2, \ldots, n} in the reference pool
#' \eqn{(y_1, \ldots, y_n)}, the predictive PMF is built from
#' \eqn{y_1, \ldots, y_{t-1}} with weight \eqn{\alpha^{t-s}} for \eqn{y_s}
#' (\eqn{s < t}).  The log-score \eqn{\log p^{<t}(y_t)} is accumulated over
#' all \eqn{t} for which \eqn{y_t} was seen at least once before \eqn{t}.
#' The total is maximised over \eqn{\alpha \in (0, 1]} via
#' \code{\link[stats]{optimize}}.
#'
#' @section \code{generate} vs \code{forecast}:
#' \describe{
#'   \item{\code{generate}}{Sample paths are drawn \emph{recursively}:
#'     \eqn{y_{T+1}^*} is sampled from the weighted pool
#'     \eqn{(y_1, \ldots, y_T)}; \eqn{y_{T+2}^*} is then sampled from the
#'     extended pool \eqn{(y_1, \ldots, y_T, y_{T+1}^*)} where the newly
#'     drawn value is the most recent and therefore receives the highest
#'     weight \eqn{\alpha^1}. Each additional horizon appends one more draw.}
#'   \item{\code{forecast}}{The predictive distribution at every horizon is
#'     the weighted empirical distribution of the reference pool with weights
#'     \eqn{\alpha^{T+h-t}} (\eqn{h}-step aging). After normalisation these
#'     weights are identical across horizons (the \eqn{\alpha^h} factor
#'     cancels), so the forecast distribution is stationary in \eqn{h}. The
#'     distributional differences across horizons are captured by the
#'     recursive \code{generate} method.}
#' }
#'
#' @param formula Model specification.
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
#' @importFrom rlang abort
#' @importFrom distributional dist_sample
#' @importFrom stats optimize
#' @export
WEMPDISTR <- function(formula, hot_start = FALSE, alpha = NULL, ...) {
  if (!is.null(alpha) &&
      (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha > 1)) {
    abort("`alpha` must be a single numeric value in (0, 1] or NULL.")
  }
  wempdistr_model <- new_model_class(
    "WEMPDISTR",
    train    = train_wempdistr,
    specials = new_specials(
      xreg = wempdistr_no_xreg
    )
  )
  new_model_definition(wempdistr_model, {{ formula }},
                       hot_start = hot_start, alpha = alpha, ...)
}


train_wempdistr <- function(.data, specials, hot_start = FALSE, alpha = NULL, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by WEMPDISTR.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by WEMPDISTR.")
  }

  T_full    <- length(y)
  start     <- if (hot_start) min(which(y != 0)) else 1L
  y_trimmed <- y[start:T_full]

  # Reference pool: last half of the (trimmed) series
  n_trimmed  <- length(y_trimmed)
  half_start <- ceiling(n_trimmed / 2)
  y_emp      <- y_trimmed[half_start:n_trimmed]

  fitted_vals <- rep(mean(y_emp), T_full)
  residuals   <- y - fitted_vals

  alpha_fit <- if (is.null(alpha)) wempdistr_loocv(y_emp) else alpha

  structure(
    list(
      y_emp     = y_emp,
      alpha     = alpha_fit,
      T_full    = T_full,
      fitted    = fitted_vals,
      residuals = residuals
    ),
    class = "WEMPDISTR"
  )
}


# Time-series leave-one-out CV for the decay parameter alpha.
#
# For each t = 2, ..., n in the reference pool y_emp, the weighted PMF is
# built from y_emp[1:(t-1)] with weight alpha^(t - s) for position s < t.
# The most-recent observation (s = t-1) therefore receives the highest weight
# alpha^1.  The log-score at y_emp[t] is accumulated; observations where
# y_emp[t] has never appeared before are skipped.
wempdistr_loocv <- function(y_emp) {
  n <- length(y_emp)
  if (n < 3L) return(0.9)   # too few points for meaningful CV

  neg_log_score <- function(alpha) {
    total <- 0
    valid <- 0L

    for (t in 2L:n) {
      hist <- y_emp[seq_len(t - 1L)]

      # raw_w[s] = alpha^(t - s) for s = 1, ..., t-1
      # s = 1 (oldest): alpha^(t-1)  |  s = t-1 (newest): alpha^1
      raw_w  <- alpha ^ seq.int(t - 1L, 1L, by = -1L)
      target <- y_emp[t]
      p      <- sum(raw_w[hist == target]) / sum(raw_w)

      if (p > 0) {
        total <- total + log(p)
        valid <- valid + 1L
      }
    }

    if (valid == 0L) return(Inf)
    -total
  }

  opt <- optimize(neg_log_score, interval = c(1e-6, 1 - 1e-6), tol = 1e-4)
  # Clamp to (0, 1] – alpha = 1 (uniform) is admissible
  min(max(opt$minimum, 1e-6), 1)
}


#' @export
forecast.WEMPDISTR <- function(object, new_data, specials = NULL,
                               n_samples = 10000L, ...) {
  h     <- nrow(new_data)
  y_emp <- object$y_emp
  alpha <- object$alpha
  n     <- length(y_emp)

  # Weights for horizon h: w_k = alpha^(n + h - k), k = 1,...,n.
  # After normalisation the alpha^h factor cancels, so the forecast
  # distribution is the same for every horizon.  We compute with the base
  # weights (h = 1) for efficiency.
  #
  # Base weights: y_emp[k] gets alpha^(n + 1 - k)
  #   k = n (y_T, most recent): alpha^1
  #   k = 1 (oldest):           alpha^n
  raw_w <- alpha ^ seq.int(n, 1L, by = -1L)

  samples <- vector("list", h)
  for (i in seq_len(h)) {
    samples[[i]] <- sample(y_emp, size = n_samples, replace = TRUE, prob = raw_w)
  }

  dist_sample(samples)
}


#' @export
generate.WEMPDISTR <- function(x, new_data, specials = NULL, ...) {
  h     <- nrow(new_data)
  alpha <- x$alpha

  # Running pool: starts as the reference pool; each draw is appended so that
  # the drawn value becomes the new most-recent observation.
  pool <- x$y_emp
  sim  <- numeric(h)

  for (i in seq_len(h)) {
    m     <- length(pool)
    # pool[k] gets weight alpha^(m + 1 - k):
    #   pool[m] (most recent): alpha^1
    #   pool[1] (oldest):      alpha^m
    raw_w  <- alpha ^ seq.int(m, 1L, by = -1L)
    sim[i] <- sample(pool, size = 1L, replace = TRUE, prob = raw_w)
    pool   <- c(pool, sim[i])
  }

  new_data$.sim <- sim
  new_data
}


#' @export
fitted.WEMPDISTR <- function(object, ...) {
  object$fitted
}


#' @export
residuals.WEMPDISTR <- function(object, ...) {
  object$residuals
}


#' @export
model_sum.WEMPDISTR <- function(x) {
  sprintf("WEMPDISTR(alpha=%.3f)", x$alpha)
}


wempdistr_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by WEMPDISTR.")
}
