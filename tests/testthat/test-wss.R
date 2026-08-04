for (i in 1:length(test_data)){
  test_that(paste0("WSS fits, forecasts, and generates on t.s. ", i), {
    test_ts <- test_data[[i]]

    # Check that the model fits correctly
    expect_no_error({
      fit <- fabletools::model(test_ts, model = WSS(value))
    })
    expect_s3_class(fit, "mdl_df")
    expect_identical(fabletools::model_sum(fit$model[[1]]), "WSS")

    # Check that fitted values and residuals are returned correctly
    fitted_vals <- stats::fitted(fit)
    resid_vals <- stats::residuals(fit)
    expect_equal(nrow(fitted_vals), nrow(test_ts))
    expect_equal(nrow(resid_vals), nrow(test_ts))

    # Check that forecasts are produced correctly
    h <- 10
    expect_no_error({
      fc <- fabletools::forecast(fit, h = h, times = 100)
    })
    expect_s3_class(fc, "fbl_ts")

    # Check the forecasts contain the expected components
    fc_mean <- fc$.mean
    fc_distr <- fc[[fabletools::distribution_var(fc)]]
    fc_family <- unname(stats::family(fc_distr))
    expect_equal(length(fc_mean), h)
    expect_equal(length(fc_distr), h)
    expect_all_true(is.finite(fc_mean))
    expect_true(inherits(fc_distr, "distribution"))
    expect_all_equal(fc_family,  "sample")
    
    # Check that simulation runs without error
    sims <- fabletools::generate(fit, h = h, times = 1)
    expect_equal(nrow(sims), h)
    expect_true(all(is.finite(sims$.sim)))

    # Check tidy
    t <- generics::tidy(fit)
    expect_s3_class(t, "tbl_df")
    expect_true(all(c("term", "estimate") %in% names(t)))
    expect_gt(nrow(t), 0L)

    # Check report
    expect_output(fabletools::report(fit))
  })
}

test_that("WSS validates inputs and rejects unsupported options", {
  expect_error(
    fable.intermittent:::train_wss(multivariate_ts(), specials = list()),
    "Only univariate responses"
  )
  expect_error(
    fable.intermittent:::train_wss(all_na_ts(), specials = list()),
    "All observations are missing"
  )
  expect_error(
    fable.intermittent:::train_wss(some_na_ts(), specials = list()),
    "Missing values are not supported"
  )
  expect_error(
    fable.intermittent:::train_wss(all_zero_ts(), specials = list()),
    "all zero"
  )
  expect_error(
    fable.intermittent:::no_xreg(),
    "Exogenous regressors are not supported"
  )

  fit <- fabletools::model(base_ts(), WSS(value))
  expect_error(
    fabletools::forecast(fit, h = 5, times = 0),
    "`times` must be a positive integer"
  )
  expect_error(
    fabletools::forecast(fit, h = 5, times = 1.5),
    "`times` must be a positive integer"
  )
})

test_that("WSS handles a series that is never zero (unvisited occurrence state)", {
  ts_never_zero <- base_ts()
  ts_never_zero$value <- abs(ts_never_zero$value) + 1
  expect_no_error({
    fit <- fabletools::model(ts_never_zero, WSS(value))
  })
  fitted_vals <- stats::fitted(fit)$.fitted
  expect_all_true(is.finite(fitted_vals[-1]))
})
