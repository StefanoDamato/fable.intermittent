for (i in 1:length(test_data)){
  test_that(paste0("GAMPOISB fits, forecasts, and generates on t.s. ", i), {
    test_ts <- test_data[[i]]
    
    # Check that the model fits correctly
    expect_no_error({
      fit <- fabletools::model(test_ts, model = GAMPOISB(value))
    })
    expect_s3_class(fit, "mdl_df")
    expect_identical(fabletools::model_sum(fit$model[[1]]), "GAMPOISB")
    
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
    expect_equal(fc_family[1], "negbin")
    expect_all_equal(fc_family[2:h],  "sample")
    
    # Check that simulation runs without error
    sims <- fabletools::generate(fit, h = h, times = 1)
    expect_equal(nrow(sims), h)
    expect_all_true(is.finite(sims$.sim))

    # Check tidy
    t <- generics::tidy(fit)
    expect_s3_class(t, "tbl_df")
    expect_true(all(c("term", "estimate") %in% names(t)))
    expect_gt(nrow(t), 0L)
  })
}
