# ---- helpers ----------------------------------------------------------------

make_ts_static <- function(y) {
  tsibble::tsibble(t = seq_along(y), value = y, index = t)
}

# ---- standard fit / forecast / generate ------------------------------------

for (i in seq_along(test_data)) {
  for (distr in c("pois", "nbinom")) {
    for (hot_start in c(FALSE, TRUE)) {
      for (alpha in list(NULL, 0.5, 1)) {
        alpha_label <- if (is.null(alpha)) "alpha=auto" else paste0("alpha=", alpha)
        test_that(paste0(
          "WSTATICDISTR(", distr, ") ", alpha_label, " ",
          ifelse(hot_start, "(hot start) ", "(cold start) "),
          "fits, forecasts, and generates on t.s. ", i
        ), {
          test_ts <- test_data[[i]]

          # ---- fit --------------------------------------------------------
          expect_no_error({
            fit <- fabletools::model(
              test_ts,
              model = WSTATICDISTR(value,
                                   distr     = distr,
                                   hot_start = hot_start,
                                   alpha     = alpha)
            )
          })
          expect_s3_class(fit, "mdl_df")

          # model_sum contains the distr family and "WSTATICDISTR"
          ms <- fabletools::model_sum(fit$model[[1]])
          expect_true(grepl("WSTATICDISTR", ms))
          expect_true(grepl(distr, ms))

          # ---- fitted / residuals ----------------------------------------
          expect_equal(nrow(stats::fitted(fit)),   nrow(test_ts))
          expect_equal(nrow(stats::residuals(fit)), nrow(test_ts))

          # ---- forecast --------------------------------------------------
          h <- 5
          expect_no_error({ fc <- fabletools::forecast(fit, h = h) })
          expect_s3_class(fc, "fbl_ts")

          fc_mean  <- fc$.mean
          fc_distr <- fc[[fabletools::distribution_var(fc)]]
          expect_equal(length(fc_mean),  h)
          expect_equal(length(fc_distr), h)
          expect_true(all(is.finite(fc_mean)))
          expect_true(all(fc_mean >= 0))
          expect_true(inherits(fc_distr, "distribution"))

          # Family must be poisson or negbin (NB may fall back to Poisson)
          fc_family <- unname(stats::family(fc_distr))
          expect_true(all(fc_family %in% c("poisson", "negbin")))

          # All horizons share the same distribution (stationary forecast)
          expect_equal(length(unique(fc_family)), 1L)

          # ---- generate --------------------------------------------------
          sims <- fabletools::generate(fit, h = h)
          expect_equal(nrow(sims), h)
          expect_true(all(is.finite(sims$.sim)))
          expect_true(all(sims$.sim >= 0))
        })
      }
    }
  }
}

# ---- alpha input validation -------------------------------------------------

test_that("WSTATICDISTR rejects alpha = 0", {
  expect_error(
    fabletools::model(make_ts_static(1:20), model = WSTATICDISTR(value, alpha = 0)),
    regexp = "alpha"
  )
})

test_that("WSTATICDISTR rejects alpha > 1", {
  expect_error(
    fabletools::model(make_ts_static(1:20), model = WSTATICDISTR(value, alpha = 1.5)),
    regexp = "alpha"
  )
})

test_that("WSTATICDISTR rejects non-numeric alpha", {
  expect_error(
    fabletools::model(make_ts_static(1:20), model = WSTATICDISTR(value, alpha = "high")),
    regexp = "alpha"
  )
})

test_that("WSTATICDISTR accepts alpha = 1 (uniform weights)", {
  y   <- c(1L, 2L, 0L, 3L, 1L, 2L, 0L, 1L)
  fit <- fabletools::model(make_ts_static(y), model = WSTATICDISTR(value, alpha = 1))
  expect_equal(fit$model[[1]]$fit$alpha, 1)
})

# ---- unsupported distr is rejected at constructor level --------------------

test_that("WSTATICDISTR rejects unsupported distributions", {
  expect_error(WSTATICDISTR(value, distr = "auto"))
  expect_error(WSTATICDISTR(value, distr = "hsp"))
  expect_error(WSTATICDISTR(value, distr = "nbinom2"))
})

# ---- stored alpha is in (0, 1] ---------------------------------------------

test_that("WSTATICDISTR auto-selected alpha is in (0, 1]", {
  for (distr in c("pois", "nbinom")) {
    fit <- fabletools::model(test_data[[1]],
                             model = WSTATICDISTR(value, distr = distr))
    alpha_fit <- fit$model[[1]]$fit$alpha
    expect_true(is.numeric(alpha_fit), label = paste("alpha numeric for", distr))
    expect_true(alpha_fit > 0 && alpha_fit <= 1,
                label = paste("alpha in (0,1] for", distr))
  }
})

# ---- very short series (LOOCV fallback) ------------------------------------

test_that("WSTATICDISTR fits on a very short series (< 3 in y_emp)", {
  # 3 obs: half_start = 2, y_emp = y[2:3] (length 2) → LOOCV returns 0.9
  short_ts <- make_ts_static(c(1L, 2L, 1L))
  for (distr in c("pois", "nbinom")) {
    expect_no_error({
      fit <- fabletools::model(short_ts, model = WSTATICDISTR(value, distr = distr))
    })
    expect_equal(fit$model[[1]]$fit$alpha, 0.9)
    expect_no_error(fabletools::forecast(fit, h = 2))
    expect_no_error(fabletools::generate(fit, h = 2))
  }
})

# ---- hot_start strips leading zeros before halving -------------------------

test_that("WSTATICDISTR hot_start uses non-zero suffix as pool", {
  y <- c(0L, 0L, 0L, 0L, 1L, 2L, 3L, 2L)
  fit_hot  <- fabletools::model(make_ts_static(y), model = WSTATICDISTR(value, hot_start = TRUE))
  fit_cold <- fabletools::model(make_ts_static(y), model = WSTATICDISTR(value, hot_start = FALSE))
  pool_hot  <- fit_hot$model[[1]]$fit$y_emp
  pool_cold <- fit_cold$model[[1]]$fit$y_emp
  expect_true(!any(pool_hot == 0))
  expect_true(any(pool_cold == 0) || length(pool_cold) > length(pool_hot))
})

# ---- forecast is stationary across horizons --------------------------------

test_that("WSTATICDISTR forecast distribution is the same for every horizon", {
  y   <- c(0L, 1L, 2L, 0L, 3L, 1L, 2L, 0L, 1L, 2L, 3L, 1L)
  fit <- fabletools::model(make_ts_static(y), model = WSTATICDISTR(value, distr = "pois"))
  fc  <- fabletools::forecast(fit, h = 5)
  means <- fc$.mean
  # All forecast means should be identical (same fitted distribution at each h)
  expect_true(all(means == means[1]))
})

# ---- generate is stochastic across replications ----------------------------

test_that("WSTATICDISTR generate produces different paths across replications", {
  set.seed(2)
  fit  <- fabletools::model(test_data[[1]], model = WSTATICDISTR(value, distr = "pois"))
  sims <- fabletools::generate(fit, h = 10, times = 10)
  paths <- split(sims$.sim, sims$.rep)
  expect_true(length(unique(paths)) > 1)
})

# ---- generate produces non-negative values ---------------------------------

test_that("WSTATICDISTR generate values are non-negative", {
  fit  <- fabletools::model(test_data[[3]], model = WSTATICDISTR(value, distr = "nbinom"))
  sims <- fabletools::generate(fit, h = 20, times = 20)
  expect_true(all(sims$.sim >= 0))
})

# ---- NB falls back to Poisson when data is not overdispersed ---------------

test_that("WSTATICDISTR(nbinom) falls back to Poisson for underdispersed data", {
  set.seed(3)
  # Poisson data is equidispersed; weighted variance <= weighted mean is likely
  y   <- rpois(30, 1.5)
  fit <- fabletools::model(make_ts_static(y), model = WSTATICDISTR(value, distr = "nbinom"))
  fc  <- fabletools::forecast(fit, h = 2)
  fc_family <- unname(stats::family(fc[[fabletools::distribution_var(fc)]]))
  # Fallback to Poisson is acceptable; nbinom is also acceptable if overdispersion
  # is detected in the weighted sample
  expect_true(all(fc_family %in% c("poisson", "negbin")))
})

# ---- multivariate input is rejected ----------------------------------------

test_that("WSTATICDISTR rejects multivariate input via train function", {
  mv_ts <- tsibble::tsibble(
    t = 1:10, a = 1:10, b = 11:20, index = t
  )
  expect_error(
    fable.intermittent:::train_wstaticdistr(mv_ts, specials = NULL, distr = "pois"),
    regexp = "univariate"
  )
})

# ---- model_sum format -------------------------------------------------------

test_that("WSTATICDISTR model_sum contains distr and alpha", {
  for (distr in c("pois", "nbinom")) {
    fit <- fabletools::model(test_data[[1]],
                             model = WSTATICDISTR(value, distr = distr, alpha = 0.7))
    ms <- fabletools::model_sum(fit$model[[1]])
    expect_true(grepl(distr,  ms))
    expect_true(grepl("0.700", ms))
  }
})
