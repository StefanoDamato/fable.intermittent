crostons_decomp <- function(y) {
  occurrence <- ifelse(y > 0, 1L, 0L)
  d_times <- which(y > 0)
  demand <- y[d_times]
  intervals <- diff(c(0, d_times))

  list(
    occurrence = occurrence,
    demand = demand,
    intervals = intervals
  )
}


get_freq <- function(.data, period = NULL, model_name = "Model") {
  period <- fabletools::get_frequencies(period, .data)
  period <- round(as.numeric(period[[1]]))

  period <- as.integer(period)
  if (period < 1) {
    rlang::abort("The seasonal period must be greater than or equal to 1.")
  }
  
  period
}


make_hurdle_shifted_distr <- function(distr, pzero){
  distr <- distributional::dist_transformed(distr, function(x) x + 1, function(x) x - 1)
  distributional::dist_inflated(distr, pzero, 0)
}



fit_nbinom <- function(y) {
  if (length(y) == 0 || all(y == 0)) {
    return(c(size = 100, prob = 1 - staticdistr_epsilon))
  }
  
  fit <- tryCatch(
    nloptr(
      x0 = c(max(mean(y), staticdistr_epsilon), 0.5),
      eval_f = function(x) -mean(dnbinom(y, x[1], x[2], log = TRUE)),
      lb = c(staticdistr_epsilon, staticdistr_epsilon),
      ub = c(Inf, 1 - staticdistr_epsilon),
      opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit) || is.null(fit$solution)) {
    mu <- mean(y)
    sigmasq <- var(y)
    if (!is.na(sigmasq) && sigmasq > mu + staticdistr_epsilon) {
      size <- (mu^2) / (sigmasq - mu)
    } else {
      size <- 100
    }
    prob <- min(size / (size + mu), 1 - staticdistr_epsilon)
    return(c(size = size, prob = prob))
  }
  
  c(size = fit$solution[1], prob = fit$solution[2])
}


