for (i in seq_along(test_data)) {
  for (distr in c("pois", "nbinom")) {
    for (hot_start in c(FALSE, TRUE)) {
      for (w in list(NULL, 0, 0.5, 1)) {
        w_label <- if (is.null(w)) "w=auto" else paste0("w=", w)
        test_that(paste0(
          "SEASSTATICDISTR(", distr, ") ", w_label, " ",
          ifelse(hot_start, "(hot start) ", "(cold start) "),
          "fits, forecasts, and generates on t.s. ", i
        ), {
          test_ts <- test_data[[i]]

          # ---- fit --------------------------------------------------------
          expect_no_error({
            fit <- fabletools::model(
              test_ts,
              model = SEASSTATICDISTR(value,
                                      distr     = distr,
                                      hot_start = hot_start,
                                      w         = w)
            )
          })
          expect_s3_class(fit, "mdl_df")

          # model_sum is a non-empty string
          ms <- fabletools::model_sum(fit$model[[1]])
          expect_true(is.character(ms) && nchar(ms) > 0)

          # ---- fitted / residuals ----------------------------------------
          fitted_vals <- stats::fitted(fit)
          resid_vals  <- stats::residuals(fit)
          expect_equal(nrow(fitted_vals), nrow(test_ts))
          expect_equal(nrow(resid_vals),  nrow(test_ts))

          # ---- forecast --------------------------------------------------
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

          # Family: "mixture" when seasonal (w > 0), plain distr otherwise
          fc_family <- unname(stats::family(fc_distr))
          # NB falls back to Poisson when data is not overdispersed, so
          # "poisson" is a valid family even when distr = "nbinom".
          if (!is.null(w) && w == 0) {
            expect_true(all(fc_family %in% c("poisson", "negbin")))
          } else {
            expect_true(all(fc_family %in% c("poisson", "negbin", "mixture")))
          }

          # ---- generate --------------------------------------------------
          sims <- fabletools::generate(fit, h = h)
          expect_equal(nrow(sims), h)
          expect_true(all(is.finite(sims$.sim)))
        })
      }
    }
  }
}

# ---- w is stored in [0, 1] ------------------------------------------------
test_that("SEASSTATICDISTR stores w in [0, 1]", {
  for (distr in c("pois", "nbinom")) {
    fit <- fabletools::model(test_data[[6]],
                             model = SEASSTATICDISTR(value, distr = distr))
    w_fit <- fit$model[[1]]$fit$w
    expect_true(is.numeric(w_fit))
    expect_true(w_fit >= 0 && w_fit <= 1)
  }
})

# ---- short series warns and falls back (w = 0, period = 1) ---------------
# A 10-obs monthly series has period = 12; 10 < 2*12=24 triggers the warning fallback.
test_that("SEASSTATICDISTR warns and falls back when n < 2*period", {
  ts_short <- tsibble::tsibble(
    time  = tsibble::yearmonth("2020 Jan") + 0:9,
    value = c(0L, 1L, 0L, 2L, 0L, 0L, 1L, 0L, 3L, 0L),
    index = "time"
  )
  for (distr in c("pois", "nbinom")) {
    expect_warning(
      fit <- fabletools::model(ts_short,
                               model = SEASSTATICDISTR(value, distr = distr)),
      regexp = "shorter than 2 \\* period"
    )
    obj <- fit$model[[1]]$fit
    expect_equal(obj$w,      0)
    expect_equal(obj$period, 1L)
  }
})

# ---- mid-range series skips LOOCV and sets w = 1/period ------------------
# A 28-obs monthly series has period = 12; 28 >= 24=2*12 but 28 < 36=3*12.
test_that("SEASSTATICDISTR skips LOOCV and sets w = 1/period when 2*period <= n < 3*period", {
  set.seed(1L)
  ts_med <- tsibble::tsibble(
    time  = tsibble::yearmonth("2020 Jan") + 0:27,
    value = c(0L, 1L, 0L, 2L, 0L, 0L, 1L, 0L, 3L, 0L, 1L, 2L,
              0L, 0L, 1L, 0L, 1L, 2L, 0L, 1L, 0L, 0L, 2L, 1L,
              1L, 0L, 2L, 0L),
    index = "time"
  )
  for (distr in c("pois", "nbinom")) {
    expect_message(
      fit <- fabletools::model(ts_med,
                               model = SEASSTATICDISTR(value, distr = distr)),
      regexp = "LOOCV is skipped"
    )
    obj <- fit$model[[1]]$fit
    expect_equal(obj$w,      1 / 12)
    expect_equal(obj$period, 12L)
  }
})

# ---- unsupported distribution is rejected at constructor level ------------
test_that("SEASSTATICDISTR rejects unsupported distributions", {
  expect_error(SEASSTATICDISTR(value, distr = "auto"))
  expect_error(SEASSTATICDISTR(value, distr = "hsp"))
  expect_error(SEASSTATICDISTR(value, distr = "nbinom2"))
})

# ---- invalid w is rejected at constructor level ---------------------------
test_that("SEASSTATICDISTR rejects invalid w", {
  expect_error(SEASSTATICDISTR(value, distr = "pois", w =  1.5), regexp = "\\[0, 1\\]")
  expect_error(SEASSTATICDISTR(value, distr = "pois", w = -0.1), regexp = "\\[0, 1\\]")
})
