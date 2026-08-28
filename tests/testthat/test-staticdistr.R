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
        expect_all_equal(fc_family, "tweedie_discrete")
      } else if (distr == "auto") {
        expect_all_true(fc_family %in% c("poisson", "negbin", "inflated",
                                         "tweedie_discrete"))
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

test_that("the Tweedie is a candidate for auto and a component of the mixture", {
  for (i in seq_along(test_data)) {
    fit <- fabletools::model(test_data[[i]], model = STATICDISTR(value, distr = "auto"))
    mdl <- fit$model[[1]]$fit
    expect_true("tweedie" %in% names(mdl$ic))
    expect_length(mdl$ic, 5)
    expect_all_true(is.finite(mdl$ic))
    expect_identical(mdl$selected_distr, names(which.min(mdl$ic)))

    fit_mix <- fabletools::model(test_data[[i]], model = STATICDISTR(value, distr = "mixture"))
    comps <- distributional::parameters(fit_mix$model[[1]]$fit$pred_distr)$dist[[1]]
    expect_length(comps, 5)
  }
})

test_that("auto ranks all five candidates on a common probability scale", {
  # The discretised Tweedie contributes a probability mass, so every candidate
  # log-likelihood is negative and the AIC values are mutually comparable.
  fit <- fabletools::model(base_ts(), model = STATICDISTR(value, distr = "auto"))
  ic <- fit$model[[1]]$fit$ic
  expect_all_true(ic > 0)
  expect_all_true(is.finite(ic))
  expect_lt(max(ic) - min(ic), 1e4)
})

test_that("the default mixture stays discrete and tweedie_discrete = FALSE relaxes it", {
  set.seed(7)
  fit <- fabletools::model(base_ts(), model = STATICDISTR(value, distr = "mixture"))
  sims <- fabletools::generate(fit, h = 200, times = 1)
  expect_all_true(sims$.sim == round(sims$.sim))

  fit_cont <- fabletools::model(
    base_ts(), model = STATICDISTR(value, distr = "mixture", tweedie_discrete = FALSE)
  )
  sims_cont <- fabletools::generate(fit_cont, h = 500, times = 1)
  expect_false(all(sims_cont$.sim == round(sims_cont$.sim)))
})

test_that("the information criteria use free-parameter counts, not parameters()", {
  # parameters() on the hurdle distributions returns the dist_inflated wrapper's
  # fields (dist, x, p), charging hsp 3 parameters when it has 2 (lambda, pzero).
  d <- fable.intermittent:::make_hurdle_shifted_distr(distributional::dist_poisson(1.2), 0.5)
  y <- c(0, 0, 3, 0, 1, 5, 0, 2)
  expect_length(distributional::parameters(d), 3)   # the upstream miscount
  expect_equal(
    fable.intermittent:::staticdistr_information(d, y, "aic", 3L) -
      fable.intermittent:::staticdistr_information(d, y, "aic", 2L),
    2
  )

  # the table covers exactly the candidates auto can fit, with the right counts
  np <- fable.intermittent:::.STATICDISTR_NPARAMS
  expect_setequal(names(np), c("pois", "hsp", "nbinom", "hsnb", "tweedie"))
  expect_identical(np[["pois"]], 1L)
  expect_identical(np[["hsp"]], 2L)
  expect_identical(np[["nbinom"]], 2L)
  expect_identical(np[["hsnb"]], 3L)
  expect_identical(np[["tweedie"]], 3L)

  # auto reports one criterion value per fitted candidate, still named
  fit <- fabletools::model(base_ts(), model = STATICDISTR(value, distr = "auto"))
  expect_setequal(names(fit$model[[1]]$fit$ic), names(np))
})

test_that("tweedie = FALSE drops the Tweedie from auto and the mixture", {
  fit <- fabletools::model(base_ts(), model = STATICDISTR(value, distr = "auto", tweedie = FALSE))
  ic <- fit$model[[1]]$fit$ic
  expect_length(ic, 4)
  expect_false("tweedie" %in% names(ic))
  expect_no_match(fabletools::model_sum(fit$model[[1]]), "tweedie")

  fit_mix <- fabletools::model(base_ts(), model = STATICDISTR(value, distr = "mixture", tweedie = FALSE))
  comps <- distributional::parameters(fit_mix$model[[1]]$fit$pred_distr)$dist[[1]]
  expect_length(comps, 4)
})

test_that("STATICDISTR validates tweedie", {
  expect_error(STATICDISTR(value, tweedie = "no"), "must be a single logical value")
  expect_error(
    STATICDISTR(value, distr = "tweedie", tweedie = FALSE),
    "incompatible with `tweedie = FALSE`"
  )
})

test_that("STATICDISTR validates tweedie_discrete", {
  expect_error(
    STATICDISTR(value, distr = "auto", tweedie_discrete = FALSE),
    "requires `tweedie_discrete = TRUE`"
  )
  expect_error(
    STATICDISTR(value, tweedie_discrete = "yes"),
    "must be a single logical value"
  )
  expect_error(
    fable.intermittent:::fit_tweedie(c(0, 1.5, 2.25), discrete = TRUE),
    "non-negative integer observations"
  )
})

test_that("the discretised Tweedie is a proper pmf and its methods agree", {
  d <- fable.intermittent:::dist_tweedie_discrete(1.5, 1.0, 1.5)
  k <- 0:80
  pk <- stats::density(d, k)[[1]]

  expect_all_true(pk >= 0)
  expect_equal(sum(pk), 1, tolerance = 1e-6)
  expect_equal(mean(d), sum(k * pk), tolerance = 1e-6)
  expect_equal(distributional::variance(d)[[1]],
               sum(k^2 * pk) - sum(k * pk)^2, tolerance = 1e-6)
  expect_equal(distributional::cdf(d, c(0, 1, 3))[[1]],
               cumsum(pk)[c(1, 2, 4)], tolerance = 1e-6)

  # density is a mass on the lattice only
  expect_equal(stats::density(d, c(-1, 0.5, 2.7))[[1]], c(0, 0, 0))

  # quantile inverts the step cdf
  probs <- c(0.1, 0.5, 0.9, 0.99)
  q <- stats::quantile(d, probs)[[1]]
  expect_all_true(q == round(q))
  expect_all_true(distributional::cdf(d, q)[[1]] >= probs - 1e-9)
  expect_all_true(distributional::cdf(d, q - 1)[[1]] < probs)

  # generate() is the rounding of the underlying Tweedie
  set.seed(42)
  g <- distributional::generate(d, 50000)[[1]]
  expect_all_true(g == round(g))
  expect_equal(mean(g), mean(d), tolerance = 0.05)
  expect_equal(mean(g == 0), pk[1], tolerance = 0.02)

  # and the log-likelihood is a probability mass: strictly negative
  expect_lt(sum(distributional::log_likelihood(d, c(0, 0, 3, 0, 1, 5, 0, 2))), 0)
})

test_that("fit_tweedie keeps the power strictly inside its bounds on count data", {
  # As the power approaches 1 the Tweedie degenerates to a Poisson, whose mass
  # lies on the integer lattice: the Lebesgue density then becomes singular and
  # the log-likelihood of integer data diverges to +Inf. The bounds prevent it.
  for (y in list(stats::rpois(60, 1.3), stats::rnbinom(60, size = 1, prob = 0.3),
                 c(rep(0, 39), 7), c(0, 0, 3, 0, 1, 5, 0, 2, 0, 0, 0, 4))) {
    params <- fable.intermittent:::fit_tweedie(y, discrete = FALSE)
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

