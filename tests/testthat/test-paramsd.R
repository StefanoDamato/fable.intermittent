for (i in 1:length(test_data)){
  for (distr in c("auto", "mixture", "pois", "nbinom", "hsp", "hsnb", "mixture")) {
    for (hot_start in c(FALSE, TRUE)) {
    test_that(paste0("PARAMSD ", "with ", distr, " distribution ",
                     ifelse(hot_start, "(hot start) ", "(cold start) "), 
                     "fits, forecasts, and generates on t.s. ", i), {
      test_ts <- test_data[[i]]
      
      # Check that the model fits correctly
      expect_no_error({
       fit <- fabletools::model(test_ts, model = PARAMSD(value,  distr = distr, hot_start = hot_start))
      })
      expect_s3_class(fit, "mdl_df")
      if (distr == "auto") {
        expect_match(fabletools::model_sum(fit$model[[1]]), "^PARAMSD\\(")
      } else {
        expect_identical(fabletools::model_sum(fit$model[[1]]), paste0("PARAMSD(", distr, ")"))
      }

      # Check that fitted values and residuals are returned correctly
      fitted_vals <- stats::fitted(fit)
      resid_vals <- stats::residuals(fit)
      expect_equal(nrow(fitted_vals), nrow(test_ts))
      expect_equal(nrow(resid_vals), nrow(test_ts))
      
      # Check that forecasts are produced correctly
      h <- 10
      expect_no_error({
       fc <- fabletools::forecast(fit, h = h)
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
      if (distr == "mixture") {
        expect_all_equal(fc_family, "mixture")
      } else if (distr == "pois") {
        expect_all_equal(fc_family, "poisson")
      } else if (distr == "nbinom") {
        expect_all_equal(fc_family, "negbin")
      } else if (distr %in% c("hsp", "hsnb")) {
        expect_all_equal(fc_family, "inflated")
      } else if (distr == "auto") {
        expect_all_true(fc_family %in% c("poisson", "negbin", "inflated"))
      }
      
      # Check that simulation runs without error
      sims <- fabletools::generate(fit, h = h)
      expect_equal(nrow(sims), h)
      expect_true(all(is.finite(sims$.sim)))

      # Check tidy
      t <- generics::tidy(fit)
      expect_s3_class(t, "tbl_df")
      expect_true(all(c("term", "estimate") %in% names(t)))

      # Check report
      expect_output(fabletools::report(fit))
      })
    }
  }
}

test_that("PARAMSD validates inputs and rejects unsupported options", {
  expect_error(
    fable.intermittent:::train_paramsd(
      multivariate_ts(), specials = list(), distr = "auto",
      hot_start = FALSE, criterion = "aic"
    ),
    "Only univariate responses"
  )
  expect_error(
    fable.intermittent:::train_paramsd(
      all_na_ts(), specials = list(), distr = "auto",
      hot_start = FALSE, criterion = "aic"
    ),
    "All observations are missing"
  )
  expect_error(
    fable.intermittent:::train_paramsd(
      some_na_ts(), specials = list(), distr = "auto",
      hot_start = FALSE, criterion = "aic"
    ),
    "Missing values are not supported"
  )
  expect_error(
    fable.intermittent:::paramsd_no_xreg(),
    "Exogenous regressors are not supported"
  )
  expect_error(
    fable.intermittent:::paramsd_information(
      distributional::dist_poisson(1), y = c(1, 2, 3), criterion = "invalid"
    ),
    "Invalid criterion"
  )
})

test_that("PARAMSD supports the bic criterion", {
  fit <- fabletools::model(base_ts(), PARAMSD(value, distr = "auto", criterion = "bic"))
  expect_s3_class(fit, "mdl_df")
  expect_output(fabletools::report(fit))
})

test_that("PARAMSD handles all-zero series (no positive demand to fit hsnb on)", {
  for (distr in c("auto", "mixture", "hsnb")) {
    fit <- fabletools::model(all_zero_ts(), PARAMSD(value, distr = distr))
    expect_s3_class(fit, "mdl_df")
    expect_all_true(is.finite(stats::fitted(fit)$.fitted))
  }
})

