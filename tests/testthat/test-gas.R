test_that("GAS supports multiple distributions", {
  bern_ts <- tsibble::tsibble(
    time = tsibble::yearmonth(as.Date("2020-01-01")) + seq_len(24) - 1,
    value = rbinom(24, size = 1, prob = 0.3),
    index = time
  )

  cases <- list(
    pois = list(data = test_data[[1]], family = "poisson"),
    nbinom = list(data = test_data[[6]], family = "negbin"),
    bern = list(data = bern_ts, family = "binomial")
  )

  for (distr in names(cases)) {
    test_ts <- cases[[distr]]$data

    fit <- fabletools::model(test_ts, model = GAS(value, distr = distr))
    expect_s3_class(fit, "mdl_df")
    expect_identical(fabletools::model_sum(fit$model[[1]]), paste0("GAS[", distr, ", RW]"))

    h <- 6
    fc <- fabletools::forecast(fit, h = h, times = 25)
    fc_distr <- fc[[fabletools::distribution_var(fc)]]
    fc_family <- unname(stats::family(fc_distr))

    expect_equal(length(fc_distr), h)
    expect_equal(fc_family[1], cases[[distr]]$family)
    expect_all_equal(fc_family[-1], "sample")

    sims <- fabletools::generate(fit, new_data = tsibble::new_data(test_ts, h))
    expect_equal(nrow(sims), h)
    expect_true(all(is.finite(sims$.sim)))
  }
})
