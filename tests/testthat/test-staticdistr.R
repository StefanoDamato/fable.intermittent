for (i in 1:length(test_data)){
  for (distr in c("auto", "mixture", "pois", "nbinom", "hsp", "hsnb", "mixture", "tweedie")) {
    for (hot_start in c(FALSE, TRUE)) {
    test_that(paste0("STATICDISTR ", "with ", distr, " distribution ",
                     ifelse(hot_start, "(hot start) ", "(cold start) "), 
                     "fits, forecasts, and generates on t.s. ", i), {
      test_ts <- test_data[[i]]
      
      # Check that the model fits correctly
      expect_no_error({
       fit <- fabletools::model(test_ts, model = STATICDISTR(value,  distr = distr, hot_start = hot_start))
      })
      expect_s3_class(fit, "mdl_df")
      if (distr == "auto") {
        expect_match(fabletools::model_sum(fit$model[[1]]), "^STATICDISTR\\(")
      } else {
        expect_identical(fabletools::model_sum(fit$model[[1]]), paste0("STATICDISTR(", distr, ")"))
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
      } else if (distr == "tweedie") {
        expect_all_equal(fc_family, "tweedie")
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

test_that("STATICDISTR validates inputs and rejects unsupported options", {
  expect_error(
    fable.intermittent:::train_staticdistr(
      multivariate_ts(), specials = list(), distr = "auto",
      hot_start = FALSE, criterion = "aic"
    ),
    "Only univariate responses"
  )
  expect_error(
    fable.intermittent:::train_staticdistr(
      all_na_ts(), specials = list(), distr = "auto",
      hot_start = FALSE, criterion = "aic"
    ),
    "All observations are missing"
  )
  expect_error(
    fable.intermittent:::train_staticdistr(
      some_na_ts(), specials = list(), distr = "auto",
      hot_start = FALSE, criterion = "aic"
    ),
    "Missing values are not supported"
  )
  expect_error(
    fable.intermittent:::staticdistr_no_xreg(),
    "Exogenous regressors are not supported"
  )
  expect_error(
    fable.intermittent:::staticdistr_information(
      distributional::dist_poisson(1), y = c(1, 2, 3), criterion = "invalid"
    ),
    "Invalid criterion"
  )
})

test_that("STATICDISTR supports the bic criterion", {
  fit <- fabletools::model(base_ts(), STATICDISTR(value, distr = "auto", criterion = "bic"))
  expect_s3_class(fit, "mdl_df")
  expect_output(fabletools::report(fit))
})

test_that("STATICDISTR never selects the Tweedie for auto or mixture", {
  for (i in seq_along(test_data)) {
    for (distr in c("auto", "mixture")) {
      fit <- fabletools::model(test_data[[i]], model = STATICDISTR(value, distr = distr))
      mdl <- fit$model[[1]]$fit
      expect_false(identical(mdl$selected_distr, "tweedie"))
      expect_false("tweedie" %in% names(mdl$ic))
      expect_no_match(fabletools::model_sum(fit$model[[1]]), "tweedie")
    }
  }
})

test_that("fit_tweedie keeps the power strictly inside its bounds on count data", {
  # As the power approaches 1 the Tweedie degenerates to a Poisson, whose mass
  # lies on the integer lattice: the Lebesgue density then becomes singular and
  # the log-likelihood of integer data diverges to +Inf. The bounds prevent it.
  for (y in list(stats::rpois(60, 1.3), stats::rnbinom(60, size = 1, prob = 0.3),
                 c(rep(0, 39), 7), c(0, 0, 3, 0, 1, 5, 0, 2, 0, 0, 0, 4))) {
    params <- fable.intermittent:::fit_tweedie(y)
    expect_gt(params[["power"]], 1.2)
    expect_lt(params[["power"]], 1.8)
    expect_gt(params[["dispersion"]], 0)
    expect_equal(params[["mean"]], mean(y))

    d <- tweedieDistr::dist_tweedie(
      params[["mean"]], params[["dispersion"]], params[["power"]]
    )
    loglik <- sum(distributional::log_likelihood(d, y))
    expect_true(is.finite(loglik))
    expect_lt(loglik, 0)
  }
})

test_that("fit_tweedie falls back gracefully on an all-zero series", {
  params <- fable.intermittent:::fit_tweedie(rep(0, 40))
  expect_gt(params[["mean"]], 0)
  expect_gt(params[["dispersion"]], 0)
  expect_gt(params[["power"]], 1)
  expect_lt(params[["power"]], 2)
})

test_that("STATICDISTR handles all-zero series (no positive demand to fit hsnb on)", {
  for (distr in c("auto", "mixture", "hsnb", "tweedie")) {
    fit <- fabletools::model(all_zero_ts(), STATICDISTR(value, distr = distr))
    expect_s3_class(fit, "mdl_df")
    expect_all_true(is.finite(stats::fitted(fit)$.fitted))
  }
})

