test_cases <- list(
    list(
      x = 0.,
      mean = 1.2,
      dispersion = 0.6,
      power = 1.8
    ),
    list(
      x = c(0.1, 1.0, 2.0, 5.0),
      mean = 3.0,
      dispersion = 4.0,
      power = 1.2
    ),
    list(
      x = c(0, 0.2, 0.8, 1.3, 2.5),
      mean = c(1.0, 1.2, 1.4, 1.6, 1.8),
      dispersion = 3.,
      power = 1.3
    ),
    list(
      x = c(0, 0.2, 1.5, 3.5),
      mean = c(1.0, 1.1, 1.4, 1.9),
      dispersion = c(0.4, 2.5, 3.5, 4.5),
      power = 1.5
    ),
    list(
      x = c(0, 0.3, 1.1, 2.4),
      mean = 0.6,
      dispersion = c(0.6, 0.9, 1.5, 4.0),
      power = c(1.2, 1.3, 1.4, 1.5)
    ),
    list(
      x = c(0.4, 1.0, 1.8, 2.6, 3.1),
      mean = c(0.1, 0.6, 1.1, 1.6, 2.1),
      dispersion = 0.6,
      power = c(1.15, 1.35, 1.55, 1.75, 1.95)
    ),
    list(
      x = c(0, 1.3, 2.8, 4.),
      mean = c(0.5, 1.5, 3.5, 5.),
      dispersion = c(0.5, 0.9, 1.6, 3.0),
      power = c(1.2, 1.5, 1.8, 1.9)
    )
  )


test_that("dist_tweedie constructs and returns core moments", {

  for (case in test_cases) {
    
    # Initialize the distribution
    y <- case$x
    mu <- case$mean
    phi <- case$dispersion
    rho <- case$power
    expect_no_error({
      distr <- dist_tweedie(mean = mu, dispersion = phi, power = rho)
    })
    expect_s3_class(distr, "distribution")
    
    # Check the moments are matched
    n <- max(length(mu), length(phi), length(rho))
    mu_long <- rep_len(mu, n)
    phi_long <- rep_len(phi, n)
    rho_long <- rep_len(rho, n)
    expect_s3_class(distr, "distribution")
    expect_equal(mean(distr), mu_long)
    expect_equal(distributional::variance(distr), phi_long * mu_long^rho_long)
    
    # Check that sampling runs without error
    samples <- generate(distr, 100)
    expect_all_equal(sapply(samples, length), 100)
    expect_all_true(sapply(samples, function(s) all(is.finite(s))))
    expect_all_true(sapply(samples, function(s) all(s >= 0)))
    
    # Check the likelihood and log-likelihood
    lik <- distributional::likelihood(distr, list(y))
    expect_length(lik, length(distr))
    expect_true(all(is.finite(lik)))
    expect_true(all(lik >= 0))
  }

})


test_that("dist_tweedie validates parameters", {
  expect_error(dist_tweedie(mean = 0), "strictly positive")
  expect_error(dist_tweedie(dispersion = 0), "strictly positive")
  expect_error(dist_tweedie(power = 1), "in \\(1, 2\\)")
})


test_that("dtweedie returns the sme results as the tweedie package", {

  for (case in test_cases) {
    y <- case$x
    mu <- case$mean
    phi <- case$dispersion
    rho <- case$power
    n <- max(length(y), length(mu), length(phi), length(rho))
    
    # Check that the function runs without error
    expect_no_error({
      density <- dtweedie(y, mean = mu, dispersion = phi, power = rho, log = FALSE)
    })
    expect_length(density, n)
    expect_true(all(is.finite(density)))
    expect_true(all(density >= 0))

    # Check the density is the same computed by the tweedie package
    if (length(case$power) == 1) {
      density_alt <- tweedie::dtweedie(y = y, xi = NULL, mu = mu, 
                                       phi = phi, power = rho, verbose = FALSE)
    } else {
      y_long <- rep_len(y, n)
      mu_long <- rep_len(mu, n)
      phi_long <- rep_len(phi, n)
      rho_long <- rep_len(rho, n)
      density_alt <- sapply(1:n, function(i) {
        tweedie::dtweedie(y = y_long[i], xi = NULL, mu = mu_long[i], 
                          phi = phi_long[i], power = rho_long[i], verbose = FALSE)
      })
    }
    expect_equal(density, density_alt)
  }
})


test_that("ptweedie returns the same results as the tweedie package", {

  for (case in test_cases) {
    q <- case$x
    mu <- case$mean
    phi <- case$dispersion
    rho <- case$power
    n <- max(length(q), length(mu), length(phi), length(rho))
    
    # Check that the function runs without error
    expect_no_error({
      cdf <- ptweedie(q, mean = mu, dispersion = phi, power = rho, log = FALSE)
    })
    expect_length(cdf, n)
    expect_true(all(is.finite(cdf)))
    expect_true(all(cdf >= 0))

    # Check the density is the same computed by the tweedie package
    if (length(case$power) == 1) {
      cdf_alt <- tweedie::ptweedie(q = q, xi = NULL, mu = mu, 
                                   phi = phi, power = rho, verbose = FALSE)
    } else {
      q_long <- rep_len(q, n)
      mu_long <- rep_len(mu, n)
      phi_long <- rep_len(phi, n)
      rho_long <- rep_len(rho, n)
      cdf_alt <- sapply(1:n, function(i) {
        tweedie::ptweedie(q = q_long[i], xi = NULL, mu = mu_long[i], 
                          phi = phi_long[i], power = rho_long[i], verbose = FALSE)
      })
    }
    expect_equal(cdf, cdf_alt)
  }
})



