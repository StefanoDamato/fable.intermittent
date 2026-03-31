for (i in 1:length(test_data)){
  for (hot_start in c(FALSE, TRUE)) {
      test_that(paste0("EMPDISTR ", ifelse(hot_start, "(hot start) ", "(cold start) "), 
                       "fits, forecasts, and generates on t.s. ", i), {
    test_ts <- test_data[[i]]
    
    # Check that the model fits correctly
    expect_no_error({
     fit <- fabletools::model(test_ts, model = EMPDISTR(value, hot_start = hot_start))
    })
    expect_s3_class(fit, "mdl_df")
    expect_identical(fabletools::model_sum(fit$model[[1]]), "EMPDISTR")
    
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
    sims <- fabletools::generate(fit, h = h)
    expect_equal(nrow(sims), h)
    expect_true(all(is.finite(sims$.sim)))
    })
  }
}
