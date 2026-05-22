for (i in seq_along(test_data)) {
  for (hot_start in c(FALSE, TRUE)) {
    for (w in list(NULL, 0, 0.5, 1)) {
      w_label <- if (is.null(w)) "w=auto" else paste0("w=", w)
      test_that(paste0(
        "SEASEMPDISTR ", w_label, " ",
        ifelse(hot_start, "(hot start) ", "(cold start) "),
        "fits, forecasts, and generates on t.s. ", i
      ), {
        test_ts <- test_data[[i]]

        # ---- fit ----------------------------------------------------------
        expect_no_error({
          fit <- fabletools::model(
            test_ts,
            model = SEASEMPDISTR(value, hot_start = hot_start, w = w)
          )
        })
        expect_s3_class(fit, "mdl_df")

        # model_sum is a non-empty string
        ms <- fabletools::model_sum(fit$model[[1]])
        expect_true(is.character(ms) && nchar(ms) > 0)

        # ---- fitted / residuals ------------------------------------------
        fitted_vals <- stats::fitted(fit)
        resid_vals  <- stats::residuals(fit)
        expect_equal(nrow(fitted_vals), nrow(test_ts))
        expect_equal(nrow(resid_vals),  nrow(test_ts))

        # ---- forecast ----------------------------------------------------
        h <- 10
        expect_no_error({
          fc <- fabletools::forecast(fit, h = h)
        })
        expect_s3_class(fc, "fbl_ts")

        fc_mean  <- fc$.mean
        fc_distr <- fc[[fabletools::distribution_var(fc)]]
        expect_equal(length(fc_mean),  h)
        expect_equal(length(fc_distr), h)
        expect_true(all(is.finite(fc_mean)))
        expect_true(inherits(fc_distr, "distribution"))
        expect_true(all(unname(stats::family(fc_distr)) == "sample"))

        # ---- generate ----------------------------------------------------
        sims <- fabletools::generate(fit, h = h)
        expect_equal(nrow(sims), h)
        expect_true(all(is.finite(sims$.sim)))
      })
    }
  }
}

# ---- w is stored in [0, 1] ------------------------------------------------
test_that("SEASEMPDISTR stores w in [0, 1]", {
  fit <- fabletools::model(test_data[[6]], model = SEASEMPDISTR(value))
  w_fit <- fit$model[[1]]$fit$w
  expect_true(is.numeric(w_fit))
  expect_true(w_fit >= 0 && w_fit <= 1)
})

# ---- short series falls back (w = 0, period = 1) -------------------------
# A 10-obs monthly series has period = 12; 10 <= 2*12 triggers the fallback.
test_that("SEASEMPDISTR falls back for short series", {
  ts_short <- tsibble::tsibble(
    time  = tsibble::yearmonth("2020 Jan") + 0:9,
    value = c(0L, 1L, 0L, 2L, 0L, 0L, 1L, 0L, 3L, 0L),
    index = "time"
  )
  fit <- fabletools::model(ts_short, model = SEASEMPDISTR(value))
  obj <- fit$model[[1]]$fit
  expect_equal(obj$w,      0)
  expect_equal(obj$period, 1L)
})

# ---- invalid w is rejected at constructor level ---------------------------
test_that("SEASEMPDISTR rejects invalid w", {
  expect_error(SEASEMPDISTR(value, w =  1.5), regexp = "\\[0, 1\\]")
  expect_error(SEASEMPDISTR(value, w = -0.1), regexp = "\\[0, 1\\]")
})
