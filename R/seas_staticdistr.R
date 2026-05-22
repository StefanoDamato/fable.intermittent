#' Seasonal Static Parametric Distribution Mixture
#'
#' A mixture of two static parametric count distributions for intermittent
#' demand forecasting. When the time series has a seasonal period `m > 1`, the
#' predictive distribution at each horizon is a weighted combination of:
#' \itemize{
#'   \item the **full distribution**: fitted via method of moments on all past
#'     observations;
#'   \item the **seasonal distribution**: fitted via method of moments on
#'     observations at the same seasonal position only.
#' }
#' The mixing weight `w` in \[0, 1\] controls the contribution of the seasonal
#' component. If `w` is not supplied it is selected automatically via
#' leave-one-out cross-validation (LOOCV) on the log-score.
#'
#' When the detected seasonal period is 1 (non-seasonal data), or when the
#' number of observations after optional leading-zero removal is not greater
#' than `2 * period`, the method falls back to a single static distribution
#' fitted on the full series (equivalent to [STATICDISTR()]).
#'
#' @section Method of moments:
#' **Poisson**: \eqn{\hat\lambda = \bar y}.
#' **Negative binomial**: \eqn{\hat{p} = \bar y / s^2},
#' \eqn{\hat{r} = \bar y^2 / (s^2 - \bar y)}, where \eqn{s^2} is the sample
#' variance. If \eqn{s^2 \le \bar y} (underdispersed), falls back to Poisson.
#'
#' @section LOOCV details:
#' Sufficient statistics (\eqn{\sum y}, \eqn{\sum y^2}) are precomputed for
#' the full series and each seasonal sub-series. Each leave-one-out
#' distribution is then refitted in O(1) via closed-form MoM updates, making
#' each optimizer evaluation O(T). When a seasonal sub-series has only one
#' observation its leave-one-out counterpart is empty; the seasonal component
#' contributes zero probability for that observation, naturally penalising
#' \eqn{w > 0}.
#'
#' @param formula Model specification.
#' @param distr Distribution family: `"pois"` (Poisson) or `"nbinom"`
#'   (negative binomial).
#' @param hot_start Logical. If `TRUE`, leading zeros are removed from the
#'   time series before fitting.
#' @param w Numeric in \[0, 1\] or `NULL` (default). The mixing weight for the
#'   seasonal distribution. When `NULL`, `w` is selected automatically via
#'   LOOCV.
#' @param ... Not used.
#'
#' @return A model specification.
#'
#' @importFrom fabletools new_model_class new_specials new_model_definition
#' @importFrom tsibble measured_vars
#' @importFrom rlang abort arg_match
#' @importFrom distributional dist_poisson dist_negative_binomial dist_mixture
#' @importFrom stats dpois dnbinom var optimize
#' @export
SEASSTATICDISTR <- function(formula,
                            distr     = c("pois", "nbinom"),
                            hot_start = FALSE,
                            w         = NULL,
                            ...) {
  distr <- arg_match(distr)
  if (!is.null(w) && (!is.numeric(w) || length(w) != 1L || w < 0 || w > 1)) {
    abort("`w` must be a single numeric value in [0, 1] or NULL.")
  }

  seasstaticdistr_model <- new_model_class(
    "SEASSTATICDISTR",
    train    = train_seasstaticdistr,
    specials = new_specials(
      xreg = seasstaticdistr_no_xreg
    )
  )
  new_model_definition(seasstaticdistr_model, {{ formula }},
                       distr = distr, hot_start = hot_start, w = w, ...)
}


train_seasstaticdistr <- function(.data, specials,
                                  distr, hot_start = FALSE, w = NULL, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by SEASSTATICDISTR.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by SEASSTATICDISTR.")
  }
  T_full <- length(y)
  start  <- if (hot_start) min(which(y != 0)) else 1L
  y_emp  <- y[start:T_full]

  full_distr  <- seasstaticdistr_mom_fit(y_emp, distr)
  fitted_vals <- rep(mean(full_distr), T_full)
  residuals   <- y - fitted_vals

  period <- get_freq(.data)

  # ---- Non-seasonal / too-short fallback ------------------------------------
  if (period <= 1L || length(y_emp) <= 2L * period) {
    return(structure(
      list(
        full_distr  = full_distr,
        seas_distrs = NULL,
        seasons_idx = NULL,
        w           = 0,
        period      = 1L,
        T_full      = T_full,
        distr       = distr,
        fitted      = fitted_vals,
        residuals   = residuals
      ),
      class = "SEASSTATICDISTR"
    ))
  }

  # ---- Season indices for y_emp ---------------------------------------------
  k_seq       <- seq_along(y_emp)
  seasons_idx <- (start + k_seq - 2L) %% period + 1L

  # ---- Seasonal distributions -----------------------------------------------
  seas_distrs <- vector("list", period)
  names(seas_distrs) <- as.character(seq_len(period))
  for (s in seq_len(period)) {
    y_s <- y_emp[seasons_idx == s]
    if (length(y_s) > 0L) {
      seas_distrs[[as.character(s)]] <- seasstaticdistr_mom_fit(y_s, distr)
    }
  }

  # ---- Mixing weight --------------------------------------------------------
  w_fit <- if (is.null(w)) {
    seasstaticdistr_loocv(y_emp, seasons_idx, period, distr)
  } else {
    w
  }

  structure(
    list(
      full_distr  = full_distr,
      seas_distrs = seas_distrs,
      seasons_idx = seasons_idx,
      w           = w_fit,
      period      = period,
      T_full      = T_full,
      distr       = distr,
      fitted      = fitted_vals,
      residuals   = residuals
    ),
    class = "SEASSTATICDISTR"
  )
}


# Fit a parametric distribution via method of moments.
seasstaticdistr_mom_fit <- function(y, distr) {
  if (distr == "pois") {
    lambda <- max(mean(y), seasstaticdistr_epsilon)
    return(dist_poisson(lambda))
  }

  if (distr == "nbinom") {
    mu <- mean(y)
    if (mu <= seasstaticdistr_epsilon || length(y) < 2L) {
      return(dist_poisson(max(mu, seasstaticdistr_epsilon)))
    }
    sigma2 <- var(y)
    if (!is.finite(sigma2) || sigma2 <= mu) {
      # Equidispersed or underdispersed: fall back to Poisson
      return(dist_poisson(mu))
    }
    size <- mu^2 / (sigma2 - mu)
    prob <- mu / sigma2
    size <- max(size, seasstaticdistr_epsilon)
    prob <- min(max(prob, seasstaticdistr_epsilon), 1 - seasstaticdistr_epsilon)
    return(dist_negative_binomial(size, prob))
  }

  abort(sprintf(
    "Distribution '%s' is not implemented in SEASSTATICDISTR. Use 'pois' or 'nbinom'.",
    distr
  ))
}


# Compute the log-density of y_t under a MoM-fitted distribution, given the
# precomputed sufficient statistics of the leave-one-out sample.
# Returns -Inf when the LOO sample is empty (n_loo <= 0).
seasstaticdistr_loo_logdens <- function(y_t, sum_loo, sumsq_loo, n_loo, distr) {
  if (n_loo <= 0L) return(-Inf)

  mu_loo <- max(sum_loo / n_loo, seasstaticdistr_epsilon)

  # Poisson, or NB with too few LOO observations to estimate variance
  if (distr == "pois" || n_loo < 2L) {
    return(dpois(y_t, mu_loo, log = TRUE))
  }

  # NB: estimate variance from LOO sufficient statistics
  # var = (sum(x^2) - n * mean^2) / (n - 1)  [Bessel-corrected]
  var_loo <- (sumsq_loo - n_loo * mu_loo^2) / (n_loo - 1L)

  if (!is.finite(var_loo) || var_loo <= mu_loo) {
    return(dpois(y_t, mu_loo, log = TRUE))
  }

  size_loo <- mu_loo^2 / (var_loo - mu_loo)
  prob_loo <- mu_loo / var_loo
  size_loo <- max(size_loo, seasstaticdistr_epsilon)
  prob_loo <- min(max(prob_loo, seasstaticdistr_epsilon), 1 - seasstaticdistr_epsilon)

  dnbinom(y_t, size_loo, prob_loo, log = TRUE)
}


# Select mixing weight w via LOOCV on the log-score.
# Sufficient statistics are precomputed so each optimizer call is O(T).
seasstaticdistr_loocv <- function(y_emp, seasons_idx, period, distr) {
  n <- length(y_emp)
  if (n < 2L) return(0)

  # Precompute sufficient statistics for the full series
  sum_y_full   <- sum(y_emp)
  sumsq_y_full <- sum(y_emp^2)

  # Precompute per-season sufficient statistics
  sum_y_seas   <- numeric(period)
  sumsq_y_seas <- numeric(period)
  n_seas_vec   <- integer(period)
  for (s in seq_len(period)) {
    idx_s          <- which(seasons_idx == s)
    n_seas_vec[s]  <- length(idx_s)
    if (length(idx_s) > 0L) {
      sum_y_seas[s]   <- sum(y_emp[idx_s])
      sumsq_y_seas[s] <- sum(y_emp[idx_s]^2)
    }
  }

  # Objective: negative total LOO log-score (minimised by optimize())
  neg_log_score <- function(w) {
    total <- 0
    valid <- 0L

    for (t in seq_len(n)) {
      y_t <- y_emp[t]
      s_t <- seasons_idx[t]

      # LOO full-distribution log-density
      ld_full <- seasstaticdistr_loo_logdens(
        y_t,
        sum_y_full - y_t,
        sumsq_y_full - y_t^2,
        n - 1L,
        distr
      )

      # LOO seasonal log-density
      # When n_s == 1, the LOO sub-series is empty → -Inf → p_seas = 0.
      # This naturally penalises w > 0 for singleton seasons.
      n_s <- n_seas_vec[s_t]
      ld_seas <- if (n_s > 1L) {
        seasstaticdistr_loo_logdens(
          y_t,
          sum_y_seas[s_t] - y_t,
          sumsq_y_seas[s_t] - y_t^2,
          n_s - 1L,
          distr
        )
      } else {
        -Inf
      }

      p_full <- exp(ld_full)
      p_seas <- if (is.finite(ld_seas)) exp(ld_seas) else 0
      p_mix  <- (1 - w) * p_full + w * p_seas

      if (is.finite(p_mix) && p_mix > 0) {
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
forecast.SEASSTATICDISTR <- function(object, new_data, specials = NULL, ...) {
  h <- nrow(new_data)

  # ---- Non-seasonal / pure-full shortcut ------------------------------------
  if (object$period <= 1L || object$w == 0) {
    return(rep(object$full_distr, h))
  }

  # Season for each forecast horizon
  next_seasons <- (object$T_full + seq_len(h) - 1L) %% object$period + 1L

  # Precompute one mixture distribution per unique season to avoid redundant
  # dist_mixture() calls when h > period.
  unique_seasons <- unique(next_seasons)
  mix_by_season  <- vector("list", object$period)

  for (s in unique_seasons) {
    seas_d <- object$seas_distrs[[as.character(s)]]
    mix_by_season[[s]] <- if (is.null(seas_d)) {
      object$full_distr
    } else {
      distributional::dist_mixture(object$full_distr, seas_d,
                                   weights = c(1 - object$w, object$w))
    }
  }

  distrs <- lapply(next_seasons, function(s) mix_by_season[[s]])
  do.call(c, distrs)
}


#' @export
generate.SEASSTATICDISTR <- function(x, new_data, specials = NULL, ...) {
  h <- nrow(new_data)

  if (x$period <= 1L || x$w == 0) {
    new_data$.sim <- unlist(distributional::generate(x$full_distr, h))
    return(new_data)
  }

  next_seasons <- (x$T_full + seq_len(h) - 1L) %% x$period + 1L
  sim <- numeric(h)

  for (i in seq_len(h)) {
    s      <- as.character(next_seasons[i])
    seas_d <- x$seas_distrs[[s]]
    d      <- if (!is.null(seas_d) && runif(1L) < x$w) seas_d else x$full_distr
    sim[i] <- distributional::generate(d, 1L)[[1L]]
  }

  new_data$.sim <- as.numeric(sim)
  new_data
}


#' @export
fitted.SEASSTATICDISTR <- function(object, ...) {
  object$fitted
}


#' @export
residuals.SEASSTATICDISTR <- function(object, ...) {
  object$residuals
}


#' @export
model_sum.SEASSTATICDISTR <- function(x) {
  if (x$period <= 1L) {
    sprintf("SEASSTATICDISTR(%s, non-seasonal)", x$distr)
  } else {
    sprintf("SEASSTATICDISTR(%s, m=%d, w=%.3f)", x$distr, x$period, x$w)
  }
}


seasstaticdistr_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by SEASSTATICDISTR.")
}

seasstaticdistr_epsilon <- 1e-4
