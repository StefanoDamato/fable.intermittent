for (i in 1:length(test_data)){
  test_that(paste0("NNARMA fits, forecasts, and generates on t.s. ", i), {
    test_ts <- test_data[[i]]
    
    # Check that the model fits correctly
    expect_no_error({
      fit <- fabletools::model(test_ts, model = NNARMA(value))
    })
    expect_s3_class(fit, "mdl_df")
    expect_identical(fabletools::model_sum(fit$model[[1]]), "NNARMA")
    
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

test_that("NNARMA validates inputs and rejects unsupported options", {
  expect_error(
    fable.intermittent:::train_nnarma(all_na_ts(), specials = list()),
    "All observations are missing"
  )
  expect_error(
    fable.intermittent:::train_nnarma(some_na_ts(), specials = list()),
    "Missing values are not supported"
  )
  expect_error(
    fable.intermittent:::train_nnarma(all_zero_ts(), specials = list()),
    "all zero"
  )
  expect_error(
    fable.intermittent:::nnarma_no_xreg(),
    "Exogenous regressors are not supported"
  )

  testthat::local_mocked_bindings(
    get_freq = function(...) 0L,
    .package = "fable.intermittent"
  )
  expect_error(
    fable.intermittent:::train_nnarma(base_ts(), specials = list()),
    "seasonal period must be greater than or equal to 1"
  )
})

test_that("NNARMA forecast floors a degenerate (zero) error variance", {
  object <- structure(
    list(
      phi = 0.5, theta = 0.1, co = 1,
      frequency = 1L, seasons = NULL,
      v_state = rep(0, 20), last_v = 0, last_m = 1, last_y = 1
    ),
    class = "NNARMA"
  )
  new_data <- tsibble::tsibble(
    time = as.Date("2026-01-01") + seq_len(5),
    index = time
  )
  fc <- fable.intermittent:::forecast.NNARMA(object, new_data)
  expect_true(all(is.finite(distributional::variance(fc))))
  expect_true(all(distributional::variance(fc) > 0))
})
