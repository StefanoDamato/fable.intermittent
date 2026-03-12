test_that("STATICDISTR fits, forecasts, and generates on all test series", {
  skip_if_not_installed("nloptr")
  set.seed(123)

  series_list <- unique(test_data$series)

  for (s in series_list) {
    data <- test_data |>
      dplyr::filter(series == s) |>
      dplyr::select(-series)

    fit <- fabletools::model(data, model = STATICDISTR(value))
    expect_s3_class(fit, "mdl_df")
    expect_identical(fabletools::model_sum(fit$model[[1]]), "STATICDISTR")

    fitted_vals <- stats::fitted(fit)
    resid_vals <- stats::residuals(fit)
    expect_equal(nrow(fitted_vals), nrow(data))
    expect_equal(nrow(resid_vals), nrow(data))

    h <- 10
    fc <- fabletools::forecast(fit, h = h, times = 100)
    expect_s3_class(fc, "fbl_ts")
    fc_mean <- fc$.mean
    fc_distr <- fc[[fabletools::distribution_var(fc)]]
    expect_equal(length(fc_mean), h)
    expect_equal(length(fc_distr), h)
    expect_true(all(is.finite(fc_mean)))
    expect_true(inherits(fc_distr, "distribution"))

    sims <- fabletools::generate(fit, h = h, times = 1)
    expect_equal(nrow(sims), h)
    expect_true(all(is.finite(sims$.sim)))
  }
})

test_that("STATICDISTR supports mixture mode", {
  data <- test_data |>
    dplyr::filter(series == "TS1") |>
    dplyr::select(-series)

  fit <- fabletools::model(data, model = STATICDISTR(value, distr = "mixture"))
  fc <- fabletools::forecast(fit, h = 6, times = 100)

  expect_equal(nrow(fc), 6)
  expect_true(inherits(fc[[fabletools::distribution_var(fc)]], "distribution"))
})
