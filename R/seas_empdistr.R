#' Seasonal Empirical Distribution Mixture
#'
#' A mixture of the full empirical distribution and a seasonal empirical
#' distribution for intermittent demand forecasting. When the time series has
#' a seasonal period `m > 1`, the predictive distribution at each horizon is a
#' weighted combination of:
#' \itemize{
#'   \item the **full empirical distribution**: built from all past observations
#'     (equivalent to [EMPDISTR()]);
#'   \item the **seasonal empirical distribution**: built only from observations
#'     at the same seasonal position (e.g., same weekday for weekly data).
#' }
#' The mixing weight `w` in \[0, 1\] controls the contribution of the seasonal
#' component (`w = 0` → pure full empirical; `w = 1` → pure seasonal). If `w`
#' is not supplied it is selected automatically via leave-one-out
#' cross-validation (LOOCV) on the log-score.
#'
#' When the detected seasonal period is 1 (non-seasonal data), or when the
#' (trimmed) series length `n` is shorter than the seasonal period, the method
#' falls back to the plain empirical distribution (a warning is issued for the
#' latter case). When `n >= period` but `n < 3 * period`, LOOCV is skipped and
#' `w` is set to `1 / period`.
#'
#' @section LOOCV details:
#' LOOCV is only run when `n >= 3 * period`. For each observation \eqn{t},
#' the leave-one-out log-score is
#' \eqn{\log\bigl[(1-w)\,p^{-t}_{\text{full}}(y_t) + w\,p^{-t}_{\text{seas}}(y_t)\bigr]},
#' where \eqn{p^{-t}} denotes the empirical probability computed without
#' observation \eqn{t}. Observations for which both leave-one-out
#' probabilities are zero (i.e., \eqn{y_t} appears only once in the full
#' series and only once in its seasonal sub-series) are skipped. The total
#' log-score is maximised over \eqn{w \in [0,1]} using [stats::optimize()].
#'
#' @param formula Model specification.
#' @param hot_start Logical. If `TRUE`, leading zeros are removed from the
#'   time series before fitting.
#' @param w Numeric in \[0, 1\] or `NULL` (default). The mixing weight for the
#'   seasonal empirical distribution. When `NULL`, `w` is selected
#'   automatically via LOOCV.
#' @param ... Not used.
#'
#' @return A model specification.
#'
#' @importFrom fabletools new_model_class new_specials new_model_definition
#' @importFrom tsibble measured_vars
#' @importFrom rlang abort
#' @importFrom distributional dist_sample
#' @importFrom stats optimize rbinom
#' @export
SEASEMPDISTR <- function(formula, hot_start = FALSE, w = NULL, ...) {
  if (!is.null(w) && (!is.numeric(w) || length(w) != 1L || w < 0 || w > 1)) {
    abort("`w` must be a single numeric value in [0, 1] or NULL.")
  }
  seasempdistr_model <- new_model_class(
    "SEASEMPDISTR",
    train   = train_seasempdistr,
    specials = new_specials(
      xreg = seasempdistr_no_xreg
    )
  )
  new_model_definition(seasempdistr_model, {{ formula }},
                       hot_start = hot_start, w = w, ...)
}


train_seasempdistr <- function(.data, specials, hot_start = FALSE, w = NULL, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by SEASEMPDISTR.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by SEASEMPDISTR.")
  }
  T_full <- length(y)
  start  <- if (hot_start) min(which(y != 0)) else 1L
  y_emp  <- y[start:T_full]

  fitted_vals <- rep(mean(y_emp), T_full)
  residuals   <- y - fitted_vals

  period <- get_freq(.data)
  n_emp  <- length(y_emp)

  # ---- Non-seasonal / short-series fallback ---------------------------------
  # Fall back when non-seasonal OR the trimmed series is shorter than one full
  # seasonal cycle (n < period): there is not enough data to build any
  # seasonal sub-series at all.
  if (period <= 1L || n_emp < period) {
    if (period > 1L) {
      warning(sprintf(
        "Series length after trimming (%d) is shorter than the seasonal period (%d). \
Falling back to non-seasonal empirical distribution.",
        n_emp, period
      ))
    }
    return(structure(
      list(
        y_emp       = y_emp,
        y_seas_list = NULL,
        seasons_idx = NULL,
        w           = 0,
        period      = 1L,
        T_full      = T_full,
        fitted      = fitted_vals,
        residuals   = residuals
      ),
      class = "SEASEMPDISTR"
    ))
  }

  # ---- Season index for every element of y_emp ------------------------------
  # y_emp[k] = y[start + k - 1]; its 1-based position in the original series
  # is (start + k - 1), giving season = (start + k - 2) %% period + 1.
  k_seq       <- seq_along(y_emp)
  seasons_idx <- (start + k_seq - 2L) %% period + 1L

  y_seas_list <- split(y_emp, seasons_idx)

  # ---- Mixing weight --------------------------------------------------------
  # Three regimes (when w is not user-supplied):
  #   n < period          -> already handled above (fallback)
  #   period <= n < 3*period -> skip LOOCV, use w = 1/period
  #   n >= 3*period          -> LOOCV
  w_fit <- if (!is.null(w)) {
    w
  } else if (n_emp < 3L * period) {
    w_default <- 1 / period
    message(sprintf(
      "Too few observations for LOOCV (n = %d < 3 * period = %d). LOOCV is skipped, set w = %.4f.",
      n_emp, 3L * period, w_default
    ))
    w_default
  } else {
    seasempdistr_loocv(y_emp, seasons_idx, period)
  }

  structure(
    list(
      y_emp       = y_emp,
      y_seas_list = y_seas_list,
      seasons_idx = seasons_idx,
      w           = w_fit,
      period      = period,
      T_full      = T_full,
      fitted      = fitted_vals,
      residuals   = residuals
    ),
    class = "SEASEMPDISTR"
  )
}


# Leave-one-out cross-validation for the mixing weight w.
# Uses precomputed count tables so that each optimizer evaluation is O(n).
seasempdistr_loocv <- function(y_emp, seasons_idx, period) {
  n <- length(y_emp)
  if (n < 2L) return(0)          # cannot do LOO with a single observation

  uniq_vals <- unique(y_emp)
  n_uniq    <- length(uniq_vals)
  val_idx   <- match(y_emp, uniq_vals)   # integer index for each observation

  # Global counts per unique value
  count_full <- tabulate(val_idx, nbins = n_uniq)

  # Per-season counts: matrix[unique_value_idx, season]
  n_seas_vec     <- integer(period)
  count_seas_mat <- matrix(0L, nrow = n_uniq, ncol = period)
  for (s in seq_len(period)) {
    idx_s          <- which(seasons_idx == s)
    n_seas_vec[s]  <- length(idx_s)
    if (length(idx_s) > 0L) {
      count_seas_mat[, s] <- tabulate(val_idx[idx_s], nbins = n_uniq)
    }
  }

  # Objective: negative total LOO log-score (to be minimised)
  neg_log_score <- function(w) {
    total <- 0
    valid <- 0L
    for (t in seq_len(n)) {
      vi  <- val_idx[t]
      s_t <- seasons_idx[t]

      # LOO full-empirical probability of y_t
      # count_full[vi] >= 1 (y_t is in y_emp), so LOO count >= 0
      p_full <- (count_full[vi] - 1L) / (n - 1L)

      # LOO seasonal probability of y_t
      # y_emp[t] is in season s_t, so count_seas_mat[vi, s_t] >= 1
      n_s    <- n_seas_vec[s_t]
      p_seas <- if (n_s > 1L) (count_seas_mat[vi, s_t] - 1L) / (n_s - 1L) else 0

      p_mix <- (1 - w) * p_full + w * p_seas
      if (p_mix > 0) {
        total <- total + log(p_mix)
        valid <- valid + 1L
      }
    }
    if (valid == 0L) return(Inf)
    -total
  }

  opt <- optimize(neg_log_score, interval = c(0, 1), tol = 1e-4)
  opt$minimum
}


#' @export
forecast.SEASEMPDISTR <- function(object, new_data, specials = NULL,
                                  n_samples = 10000L, ...) {
  h <- nrow(new_data)

  # ---- Non-seasonal / pure-full-empirical shortcut --------------------------
  if (object$period <= 1L || object$w == 0) {
    return(dist_sample(rep(list(object$y_emp), h)))
  }

  # Season for each forecast horizon.
  # Forecast step j (j = 1,...,h) corresponds to 1-based position T_full + j
  # in the original series, giving season = (T_full + j - 1) %% period + 1.
  next_seasons <- (object$T_full + seq_len(h) - 1L) %% object$period + 1L
  w            <- object$w

  samples <- vector("list", h)

  for (i in seq_len(h)) {
    s   <- as.character(next_seasons[i])
    y_s <- object$y_seas_list[[s]]

    if (is.null(y_s) || length(y_s) == 0L) {
      # Season never observed in training: fall back to full empirical
      samples[[i]] <- object$y_emp
    } else {
      # Deterministic split to keep forecasts reproducible
      n_from_seas <- round(w * n_samples)
      n_from_full <- n_samples - n_from_seas

      samp_seas <- if (n_from_seas > 0L) sample(y_s,          n_from_seas, replace = TRUE) else numeric(0L)
      samp_full <- if (n_from_full > 0L) sample(object$y_emp, n_from_full, replace = TRUE) else numeric(0L)
      samples[[i]] <- c(samp_seas, samp_full)
    }
  }

  dist_sample(samples)
}


#' @export
generate.SEASEMPDISTR <- function(x, new_data, specials = NULL, ...) {
  h <- nrow(new_data)

  if (x$period <= 1L || x$w == 0) {
    sim <- sample(x$y_emp, size = h, replace = TRUE)
  } else {
    next_seasons <- (x$T_full + seq_len(h) - 1L) %% x$period + 1L
    sim          <- numeric(h)
    for (i in seq_len(h)) {
      s   <- as.character(next_seasons[i])
      y_s <- x$y_seas_list[[s]]
      if (is.null(y_s) || length(y_s) == 0L) {
        sim[i] <- sample(x$y_emp, 1L, replace = TRUE)
      } else {
        pool   <- if (runif(1L) < x$w) y_s else x$y_emp
        sim[i] <- sample(pool, 1L, replace = TRUE)
      }
    }
  }

  new_data$.sim <- as.numeric(sim)
  new_data
}


#' @export
fitted.SEASEMPDISTR <- function(object, ...) {
  object$fitted
}


#' @export
residuals.SEASEMPDISTR <- function(object, ...) {
  object$residuals
}


#' @export
model_sum.SEASEMPDISTR <- function(x) {
  if (x$period <= 1L) {
    "SEASEMPDISTR(non-seasonal)"
  } else {
    sprintf("SEASEMPDISTR(m=%d, w=%.3f)", x$period, x$w)
  }
}


seasempdistr_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by SEASEMPDISTR.")
}
