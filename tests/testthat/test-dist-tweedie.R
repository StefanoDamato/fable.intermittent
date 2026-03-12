test_that("dist_tweedie constructs and returns core moments", {
  d <- dist_tweedie(mean = 2, dispersion = 0.8, power = 1.5)

  expect_s3_class(d, "distribution")
  expect_equal(distributional::covariance(d), 0.8 * (2^1.5))
})

test_that("dist_tweedie generates finite non-negative samples", {
  set.seed(123)
  d <- dist_tweedie(mean = 1.5, dispersion = 0.6, power = 1.4)

  x <- generate(d, 200)[[1]]
  expect_equal(length(x), 200)
  expect_true(all(is.finite(x)))
  expect_true(all(x >= 0))
})

test_that("dist_tweedie validates parameters", {
  expect_error(dist_tweedie(mean = 0), "strictly positive")
  expect_error(dist_tweedie(dispersion = 0), "strictly positive")
  expect_error(dist_tweedie(power = 1), "in \\(1, 2\\)")
})
