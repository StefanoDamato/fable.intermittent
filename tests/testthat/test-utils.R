test_that("get_freq errors when the seasonal period is below 1", {
  ts <- base_ts()
  expect_error(
    fable.intermittent:::get_freq(ts, period = 0),
    "seasonal period must be greater than or equal to 1"
  )
})

test_that("dist_normal_nonneg errors on non-positive sigma", {
  expect_error(
    fable.intermittent:::dist_normal_nonneg(mu = 0, sigma = 0),
    "sigma parameter"
  )
  expect_error(
    fable.intermittent:::dist_normal_nonneg(mu = 0, sigma = -1),
    "sigma parameter"
  )
})

test_that("dist_normal_nonneg S3 methods behave as documented", {
  d <- fable.intermittent:::dist_normal_nonneg(mu = 2, sigma = 1.5)

  expect_match(format(d), "^N\\+\\(")

  dens <- unlist(stats::density(d, at = c(-1, 0, 1)))
  expect_equal(dens[1], 0)
  expect_equal(dens[2], stats::pnorm(0, 2, 1.5))
  expect_equal(dens[3], stats::dnorm(1, 2, 1.5))

  cdf_default <- unlist(distributional::cdf(d, q = 1))
  cdf_upper <- unlist(distributional::cdf(d, q = 1, lower.tail = FALSE))
  cdf_log <- unlist(distributional::cdf(d, q = 1, log.p = TRUE))
  expect_equal(cdf_default + cdf_upper, 1)
  expect_equal(cdf_log, log(cdf_default))

  q_default <- unlist(stats::quantile(d, p = 0.9))
  q_upper <- unlist(stats::quantile(d, p = 0.1, lower.tail = FALSE))
  q_log <- unlist(stats::quantile(d, p = log(0.9), log.p = TRUE))
  expect_equal(q_default, q_upper)
  expect_equal(q_default, q_log)

  expect_equal(unlist(distributional::covariance(d)), 1.5^2)
  expect_equal(mean(d), 2)
})

test_that("fit_nbinom falls back to a moment-based estimate when optimisation fails", {
  testthat::local_mocked_bindings(
    nloptr = function(...) stop("forced optimisation failure"),
    .package = "fable.intermittent"
  )

  # Overdispersed data: falls back to the method-of-moments size estimate
  y_overdispersed <- stats::rnbinom(30, size = 1, prob = 0.3)
  res_over <- fable.intermittent:::fit_nbinom(y_overdispersed)
  expect_named(res_over, c("size", "prob"))
  expect_all_true(is.finite(res_over))

  # Equidispersed data (variance <= mean): falls back to a fixed large size
  y_equidispersed <- rep(5, 30)
  res_equi <- fable.intermittent:::fit_nbinom(y_equidispersed)
  expect_named(res_equi, c("size", "prob"))
  expect_equal(unname(res_equi["size"]), 100)
})
