for (i in 1:length(test_data)){
  test_that(paste0("MARWAL fits, forecasts, and generates on t.s. ", i), {
    test_ts <- test_data[[i]]
    
    # Check that the model fits correctly
    expect_no_error({
     fit <- fabletools::model(test_ts, model = MARWAL(value))
    })
    expect_s3_class(fit, "mdl_df")
    expect_identical(fabletools::model_sum(fit$model[[1]]), "MARWAL")
    
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
    expect_all_equal(fc_family, "normal_nonneg")
    
    # Check that simulation runs without error
    sims <- fabletools::generate(fit, h = h)
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

test_that("MARWAL validates inputs and rejects unsupported options", {
  expect_error(
    fable.intermittent:::train_marwal(all_na_ts(), specials = list()),
    "All observations are missing"
  )
  expect_error(
    fable.intermittent:::train_marwal(some_na_ts(), specials = list()),
    "Missing values are not supported"
  )
  expect_error(
    fable.intermittent:::train_marwal(all_zero_ts(), specials = list()),
    "all zero"
  )
  expect_error(
    fable.intermittent:::marwal_no_xreg(),
    "Exogenous regressors are not supported"
  )

  testthat::local_mocked_bindings(
    get_freq = function(...) 0L,
    .package = "fable.intermittent"
  )
  expect_error(
    fable.intermittent:::train_marwal(base_ts(), specials = list()),
    "seasonal period must be greater than or equal to 1"
  )
})

test_that("MARWAL handles Markov chains with an unvisited occurrence state", {
  # Never zero: the "previous state = 0" row of the transition count is
  # never observed, exercising the `lambda <- 1` fallback.
  ts_never_zero <- base_ts()
  ts_never_zero$value <- abs(ts_never_zero$value) + 1
  expect_no_error({
    fit <- fabletools::model(ts_never_zero, MARWAL(value))
  })
  expect_all_true(is.finite(stats::fitted(fit)$.fitted))

  # Zero everywhere except the very last observation: the "previous state
  # = 1" row is never observed, exercising the `lambda <- .MARWAL_EPSILON`
  # fallback.
  ts_zero_then_one <- base_ts()
  ts_zero_then_one$value <- 0
  ts_zero_then_one$value[nrow(ts_zero_then_one)] <- 5
  expect_no_error({
    fit <- fabletools::model(ts_zero_then_one, MARWAL(value))
  })
  expect_all_true(is.finite(stats::fitted(fit)$.fitted))
})
