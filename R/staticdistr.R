#' Static Count Distribution Model
#'
#' Static (IID) count distribution model for intermittent demand, following
#' Kolassa (2016). The method fits several candidate distributions --- Poisson,
#' hurdle-shifted Poisson, negative binomial, and hurdle-shifted negative
#' binomial --- to the observed series and selects the best by AIC. A mixture
#' option that blends all four predictive distributions is also available.
#'
#' A static Tweedie distribution (`distr = "tweedie"`) is also available, and
#' takes part in `"auto"` and `"mixture"` as a fifth candidate. A Tweedie is
#' continuous on the positive half-line with an atom at zero, so its
#' log-likelihood is a density rather than a probability mass and cannot be
#' ranked against the count candidates by AIC/BIC, nor coherently blended with
#' them. `tweedie_discrete` resolves this by rounding the Tweedie to the
#' non-negative integers, giving the proper probability mass function
#' `P(Y = 0) = F(0.5)` and `P(Y = k) = F(k + 0.5) - F(k - 0.5)`.
#'
#' `distr = "auto"` therefore always uses the discretised form, and rejects
#' `tweedie_discrete = FALSE`. `distr = "mixture"` and `distr = "tweedie"`
#' honour the argument, so the continuous Tweedie remains reachable for
#' intermittent series that are not counts --- at the cost, for the mixture, of
#' blending a density with four probability masses.
#'
#' @param formula Model specification.
#' @param distr Distribution choice: one of `"auto"`, `"pois"`, `"hsp"`,
#'   `"nbinom"`, `"hsnb"`, `"mixture"`, or `"tweedie"`.
#' @param tweedie Logical. If `TRUE` (the default) the Tweedie takes part in
#'   `distr = "auto"` and `distr = "mixture"` as a fifth candidate. Set it to
#'   `FALSE` to rank and blend only the four count distributions.
#' @param tweedie_discrete Logical. If `TRUE` (the default) the Tweedie
#'   candidate is discretised by rounding to the non-negative integers, keeping
#'   it on the same probability scale as the count distributions. Ignored unless
#'   the Tweedie is fitted, and required to be `TRUE` when `distr = "auto"`.
#' @param hot_start Logical. If `TRUE`, leading zeros are removed from the
#'   time series before fitting.
#' @param criterion Information criterion to use for model selection when `distr =
#'   "auto"`. One of `"aic"` or `"bic"`.
#' @param object A fitted model object.
#' @param ... Not used.
#'
#' @references
#' Kolassa, S. (2016). Evaluating predictive count data distributions in retail
#' sales forecasting. *International Journal of Forecasting*, 32(3), 788--803.
#' \doi{10.1016/j.ijforecast.2015.12.004}.
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
#'   model(STATICDISTR(value)) |>
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
#' @importFrom rlang abort arg_match is_integerish
#' @importFrom distributional dist_poisson dist_negative_binomial log_likelihood parameters dist_sample
#' @importFrom nloptr nloptr
#' @importFrom stats dpois dnbinom rpois rnbinom runif var setNames
#' @importFrom tweedieDistr dist_tweedie
#' @export
STATICDISTR <- function(formula, distr = c("auto", "pois", "hsp", "nbinom", "hsnb",
                                           "mixture", "tweedie"),
                        hot_start = FALSE, criterion = c("aic", "bic"),
                        tweedie = TRUE, tweedie_discrete = TRUE, ...) {
  distr <- arg_match(distr)
  criterion <- arg_match(criterion)

  if (!is.logical(tweedie) || length(tweedie) != 1L || is.na(tweedie)) {
    abort("`tweedie` must be a single logical value.")
  }
  if (!is.logical(tweedie_discrete) || length(tweedie_discrete) != 1L ||
      is.na(tweedie_discrete)) {
    abort("`tweedie_discrete` must be a single logical value.")
  }
  if (distr == "tweedie" && !tweedie) {
    abort("`distr = \"tweedie\"` is incompatible with `tweedie = FALSE`.")
  }
  if (distr == "auto" && !tweedie_discrete) {
    abort(paste0(
      "`distr = \"auto\"` requires `tweedie_discrete = TRUE`: the continuous ",
      "Tweedie log-likelihood is a density, so it cannot be ranked against the ",
      "count distributions by AIC/BIC."
    ))
  }

  staticdistr_model <- new_model_class(
    "STATICDISTR",
    train = train_staticdistr,
    specials = new_specials(
      xreg = staticdistr_no_xreg
    )
  )
  new_model_definition(staticdistr_model, {{ formula }}, distr = distr,
                       hot_start = hot_start, criterion = criterion,
                       tweedie = tweedie, tweedie_discrete = tweedie_discrete, ...)
}

train_staticdistr <- function(.data, specials, distr, hot_start, criterion,
                              tweedie = TRUE, tweedie_discrete = TRUE, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by STATICDISTR.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  if (anyNA(y)) {
    abort("Missing values are not supported by STATICDISTR.")
  }

  if (hot_start) {
    start <- which(y > 0)[1]
    y <- y[start:length(y)]
  } else {
    start <- 1
  }

  # Identify the distributions to be fitted. The Tweedie joins the four count
  # candidates; whether it enters in its discretised form is decided below.
  if (distr %in% c("auto", "mixture")) {
    to_eval <- c("nbinom", "pois", "hsnb", "hsp")
    if (tweedie) {
      to_eval <- c(to_eval, "tweedie")
    }
  } else {
    to_eval <- distr
  }

  # Ranking by AIC/BIC is only meaningful on a common probability scale, so the
  # discretised Tweedie is mandatory for "auto". The mixture and the explicit
  # choice leave it to the caller.
  use_discrete <- if (distr == "auto") TRUE else tweedie_discrete


  # Apply Croston's decomposition
  decomp <- crostons_decomp(y)
  occurrence <- decomp$occurrence
  shifted_demand <- decomp$demand - 1

  #Fit the distributions
  fit_distr <- list()
  if ("pois" %in% to_eval) {
    fit_distr[["pois"]] <- staticdistr_fit_pois(y)
  }
  if ("hsp" %in% to_eval) {
    fit_distr[["hsp"]] <- staticdistr_fit_hsp(occurrence, shifted_demand)
  }
  if ("nbinom" %in% to_eval) {
    fit_distr[["nbinom"]] <- staticdistr_fit_nbinom(y)
  }
  if ("hsnb" %in% to_eval) {
    fit_distr[["hsnb"]] <- staticdistr_fit_hsnb(occurrence, shifted_demand)
  }
  if ("tweedie" %in% to_eval) {
    fit_distr[["tweedie"]] <- staticdistr_fit_tweedie(y, discrete = use_discrete)
  }

  # Select the distribution to use for forecasting
  if (distr == "mixture") {
    w <- rep(1/length(fit_distr), length(fit_distr))
    pred_distr <- do.call(distributional::dist_mixture, c(fit_distr, list(weights = w)))
    ic <- NULL
  } else if (distr == "auto") {
    ic <- vapply(names(fit_distr), function(nm) {
      staticdistr_information(fit_distr[[nm]], y, criterion, .STATICDISTR_NPARAMS[[nm]])
    }, numeric(1))
    pred_distr <- fit_distr[[names(which.min(ic))]]
  } else {
    pred_distr <- fit_distr[[distr]]
    ic <- NULL
  }
  selected_distr <- if (distr == "auto") names(which.min(ic)) else distr

  # Compute fitted values and residuals
  init_na <- rep(NA, start - 1)
  fitted <- c(init_na, rep(mean(pred_distr), length(y)))
  residuals <- c(init_na, y) - fitted

  structure(
    list(
      selected_distr = selected_distr,
      ic = ic,
      pred_distr = pred_distr,
      fitted = fitted,
      residuals = residuals
    ),
    class = "STATICDISTR"
  )
}

#' Forecast a STATICDISTR model
#'
#' Produces forecast distributions from a fitted STATICDISTR model.
#'
#' @inheritParams forecast.EMPDISTR
#'
#' @return A distribution vector. The class depends on the static distribution
#'    fitted by the `STATICDISTR` method.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, STATICDISTR(value))
#' forecast(fit, h = "7 days")
#'
#' @export
forecast.STATICDISTR <- function(object, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  rep(object$pred_distr, h)
}

#' Generate sample paths from a STATICDISTR model
#'
#' @param x A fitted `STATICDISTR` model object.
#' @inheritParams forecast.STATICDISTR
#'
#' @return A vector of future paths from a dataset using a fitted model.
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, STATICDISTR(value))
#' generate(fit, new_data = tsibble::new_data(ts, 7))
#' @export
generate.STATICDISTR <- function(x, new_data, specials = NULL, ...) {
  h <- nrow(new_data)
  new_data$.sim <- unlist(distributional::generate(x$pred_distr, h))
  new_data
}

#' Extract fitted values from a STATICDISTR model
#'
#' @inherit fitted.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, STATICDISTR(value))
#' fitted(fit)
#' @export
fitted.STATICDISTR <- function(object, ...) {
  object$fitted
}

#' Extract residuals from a STATICDISTR model
#'
#' @inherit residuals.EMPDISTR
#'
#' @examples
#' ts <- tsibble::tsibble(
#'   time = as.Date("2026-01-01") + seq_len(40),
#'   value = rnbinom(40, size = 1, prob = 0.3),
#'   index = time
#' )
#' fit <- model(ts, STATICDISTR(value))
#' residuals(fit)
#' @export
residuals.STATICDISTR <- function(object, ...) {
  object$residuals
}


#' @export
model_sum.STATICDISTR <- function(x) {
  paste0("STATICDISTR(", x$selected_distr, ")")
}

#' @export
tidy.STATICDISTR <- function(x, ...) {
  tryCatch({
    params <- as.list(distributional::parameters(x$pred_distr))
    tibble(
      term     = names(params),
      estimate = vapply(params, function(v) {
        v <- unlist(v)
        if (is.numeric(v) && length(v) == 1L) v else NA_real_
      }, numeric(1))
    )
  }, error = function(e) tibble(term = character(), estimate = numeric()))
}

#' @rdname STATICDISTR
#' @export
report.STATICDISTR <- function(object, ...) {
  tryCatch({
    params <- as.list(distributional::parameters(object$pred_distr))
    if (length(params) > 0) {
      cat("  Parameters:\n")
      for (nm in names(params)) {
        val <- unlist(params[[nm]])
        if (is.numeric(val) && length(val) == 1)
          cat(sprintf("    %-8s = %g\n", nm, val))
      }
    }
  }, error = function(e) NULL)
  if (!is.null(object$ic)) {
    cat("\n  Information criteria:\n")
    for (nm in names(object$ic))
      cat(sprintf("    %-8s = %.2f\n", nm, object$ic[[nm]]))
  }
  invisible(object)
}


staticdistr_fit_pois <- function(y) {
  lambda <- mean(y)
  distributional::dist_poisson(lambda)
}

staticdistr_fit_hsp <- function(occurrence, shifted_demand) {
  pzero = mean(1 - occurrence)
  lambda = ifelse(length(shifted_demand) > 0, mean(shifted_demand), 0)
  make_hurdle_shifted_distr(dist_poisson(lambda), pzero)
}

staticdistr_fit_nbinom <- function(y) {
  params <- fit_nbinom(y)
  distributional::dist_negative_binomial(params[['size']], params[['prob']])
}

staticdistr_fit_hsnb <- function(occurrence, shifted_demand) {
  if (length(shifted_demand) > 0) {
    params <- fit_nbinom(shifted_demand)
  } else {
    params <- c(size = 100, prob = 1 - .STATICDISTR_EPSILON)
  }
  pzero = mean(1 - occurrence)
  make_hurdle_shifted_distr(dist_negative_binomial(params[['size']], params[['prob']]), pzero)
}

staticdistr_fit_tweedie <- function(y, discrete = TRUE) {
  params <- fit_tweedie(y, discrete = discrete)
  if (discrete) {
    dist_tweedie_discrete(params[['mean']], params[['dispersion']], params[['power']])
  } else {
    dist_tweedie(params[['mean']], params[['dispersion']], params[['power']])
  }
}


# Free parameters per candidate. Deliberately not derived from parameters():
# for the hurdle distributions that returns the dist_inflated wrapper's fields
# (dist, x, p), which counts the fixed inflation point x = 0 and misses the
# inner distribution's parameters -- so hsp is charged 3 instead of 2. See
# https://github.com/mitchelloharawild/distributional/issues/161
.STATICDISTR_NPARAMS <- c(pois = 1L, hsp = 2L, nbinom = 2L, hsnb = 3L, tweedie = 3L)

# n_params defaults to the (unreliable) introspection so that direct calls
# without a candidate name keep working.
staticdistr_information <- function(distr, y, criterion, n_params = NULL){
  loglik <- sum(distributional::log_likelihood(distr, y))
  n_obs <- length(y)
  if (is.null(n_params)) {
    n_params <- length(distributional::parameters(distr))
  }

  if (criterion == "aic") {
    -2 * loglik + 2 * n_params
  } else if (criterion == "bic") {
    -2 * loglik + log(n_obs) * n_params
  } else {
    abort("Invalid criterion. Use 'aic' or 'bic'.")
  }
}


staticdistr_no_xreg <- function(...) {
  abort("Exogenous regressors are not supported by STATICDISTR.")
}

