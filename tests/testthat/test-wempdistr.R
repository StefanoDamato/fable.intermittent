# ---- helpers ----------------------------------------------------------------

make_ts <- function(y) {
  tsibble::tsibble(t = seq_along(y), value = y, index = t)
}

# ---- standard fit / forecast / generate ------------------------------------

for (i in seq_along(test_data)) {
  for (hot_start in c(FALSE, TRUE)) {
    for (alpha in list(NULL, 0.5, 1)) {
      alpha_label <- if (is.null(alpha)) "alpha=auto" else paste0("alpha=", alpha)
      test_that(paste0(
        "WEMPDISTR ", alpha_label, " ",
        ifelse(hot_start, "(hot start) ", "(cold start) "),
        "fits, forecasts, and generates on t.s. ", i
      ), {
        test_ts <- test_data[[i]]

        expect_no_error({
          fit <- fabletools::model(
            test_ts,
            model = WEMPDISTR(value, hot_start = hot_start, alpha = alpha)
          )
        })
        expect_s3_class(fit, "mdl_df")

        # model_sum contains "WEMPDISTR"
        ms <- fabletools::model_sum(fit$model[[1]])
        expect_true(grepl("WEMPDISTR", ms))

        # fitted / residuals have one row per observation
        expect_equal(nrow(stats::fitted(fit)),   nrow(test_ts))
        expect_equal(nrow(stats::residuals(fit)), nrow(test_ts))

        # forecast
        h <- 5
        expect_no_error({ fc <- fabletools::forecast(fit, h = h) })
        expect_s3_class(fc, "fbl_ts")
        expect_equal(length(fc$.mean), h)
        expect_true(all(is.finite(fc$.mean)))
        fc_distr <- fc[[fabletools::distribution_var(fc)]]
        expect_true(all(unname(stats::family(fc_distr)) == "sample"))

        # generate
        expect_no_error({ sims <- fabletools::generate(fit, h = h) })
        expect_equal(nrow(sims), h)
        expect_true(all(is.finite(sims$.sim)))
      })
    }
  }
}

# ---- alpha input validation -------------------------------------------------

test_that("WEMPDISTR rejects alpha = 0", {
  expect_error(
    fabletools::model(make_ts(1:20), model = WEMPDISTR(value, alpha = 0)),
    regexp = "alpha"
  )
})

test_that("WEMPDISTR rejects alpha > 1", {
  expect_error(
    fabletools::model(make_ts(1:20), model = WEMPDISTR(value, alpha = 1.5)),
    regexp = "alpha"
  )
})

test_that("WEMPDISTR rejects non-numeric alpha", {
  expect_error(
    fabletools::model(make_ts(1:20), model = WEMPDISTR(value, alpha = "high")),
    regexp = "alpha"
  )
})

test_that("WEMPDISTR accepts alpha = 1 (uniform weights)", {
  fit <- fabletools::model(make_ts(c(0L, 1L, 2L, 1L, 0L, 3L, 1L, 2L)), model = WEMPDISTR(value, alpha = 1))
  expect_equal(fit$model[[1]]$fit$alpha, 1)
})

# ---- stored alpha is in (0, 1] ----------------------------------------------

test_that("WEMPDISTR auto-selected alpha is in (0, 1]", {
  fit <- fabletools::model(test_data[[1]], model = WEMPDISTR(value))
  alpha_fit <- fit$model[[1]]$fit$alpha
  expect_true(is.numeric(alpha_fit))
  expect_true(alpha_fit > 0 && alpha_fit <= 1)
})

# ---- very short series (LOOCV fallback) -------------------------------------

test_that("WEMPDISTR fits on a very short series (< 3 in y_emp)", {
  # 3 observations: half_start = ceiling(3/2) = 2, so y_emp = y[2:3] (length 2)
  # → LOOCV n < 3 guard fires and returns the default 0.9
  short_ts <- make_ts(c(0L, 1L, 2L))
  expect_no_error({
    fit <- fabletools::model(short_ts, model = WEMPDISTR(value))
  })
  expect_equal(fit$model[[1]]$fit$alpha, 0.9)
  expect_no_error(fabletools::forecast(fit, h = 2))
  expect_no_error(fabletools::generate(fit, h = 2))
})

# ---- hot_start strips leading zeros before halving --------------------------

test_that("WEMPDISTR hot_start uses non-zero suffix for the pool", {
  # 8 obs: 4 leading zeros then 1,2,3,2 → after hot_start y_trimmed = 1,2,3,2
  # → y_emp = last half = 3,2  (not the zeros)
  y <- c(0L, 0L, 0L, 0L, 1L, 2L, 3L, 2L)
  fit_hot  <- fabletools::model(make_ts(y), model = WEMPDISTR(value, hot_start = TRUE))
  fit_cold <- fabletools::model(make_ts(y), model = WEMPDISTR(value, hot_start = FALSE))
  pool_hot  <- fit_hot$model[[1]]$fit$y_emp
  pool_cold <- fit_cold$model[[1]]$fit$y_emp
  # hot_start pool must contain no zeros (they were trimmed before halving)
  expect_true(!any(pool_hot == 0))
  # cold_start pool may contain zeros
  expect_true(any(pool_cold == 0) || length(pool_cold) > length(pool_hot))
})

# ---- generate is stochastic across replications -----------------------------

test_that("WEMPDISTR generate produces different paths across replications", {
  set.seed(1)
  fit  <- fabletools::model(test_data[[1]], model = WEMPDISTR(value))
  sims <- fabletools::generate(fit, h = 10, times = 10)
  paths <- split(sims$.sim, sims$.rep)
  # With high probability, not all 10 paths are identical
  expect_true(length(unique(paths)) > 1)
})

# ---- generate draws only from the support of y_emp (possibly extended) -----

test_that("WEMPDISTR generate values are within training-data support", {
  y    <- c(0L, 0L, 1L, 2L, 0L, 3L, 1L, 2L, 0L, 1L)
  fit  <- fabletools::model(make_ts(y), model = WEMPDISTR(value))
  sims <- fabletools::generate(fit, h = 20, times = 20)
  pool <- fit$model[[1]]$fit$y_emp
  expect_true(all(sims$.sim %in% pool))
})

# ---- multivariate input is rejected ----------------------------------------
# fabletools wraps training errors as warnings, so we call the internal
# training function directly to verify the guard fires.

test_that("WEMPDISTR rejects multivariate input via train function", {
  mv_ts <- tsibble::tsibble(
    t = 1:10, a = 1:10, b = 11:20, index = t
  )
  expect_error(
    fable.intermittent:::train_wempdistr(mv_ts, specials = NULL),
    regexp = "univariate"
  )
})

# ---- forecast dist family is always "sample" --------------------------------

test_that("WEMPDISTR forecast distribution family is 'sample'", {
  fit <- fabletools::model(test_data[[3]], model = WEMPDISTR(value))
  fc  <- fabletools::forecast(fit, h = 6)
  distr <- fc[[fabletools::distribution_var(fc)]]
  expect_true(all(unname(stats::family(distr)) == "sample"))
})
