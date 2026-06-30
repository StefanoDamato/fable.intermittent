################################################################################
# GLOBAL PARAMETERS (EPSILON) TO AVOID NUMERICAL ISSUES IN COMPUTATIONS
#' @importFrom distributional cdf dist_inflated dist_transformed dist_truncated
#' @importFrom fabletools get_frequencies
#' @importFrom rlang abort
NULL

.BETANBB_EPSILON     <- 1e-4
.GAMPOISB_EPSILON    <- 1e-4
.HSPES_EPSILON       <- 1e-4
.MARWAL_EPSILON      <- 1e-4
.NNARMA_EPSILON      <- 1e-4
.NEGBINES_EPSILON    <- 1e-4
.STATICDISTR_EPSILON <- 1e-4
.TWEES_EPSILON       <- 1e-4

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
  period <- get_frequencies(period, .data)
  period <- round(as.numeric(period[[1]]))

  period <- as.integer(period)
  if (period < 1) {
    abort("The seasonal period must be greater than or equal to 1.")
  }

  period
}

make_hurdle_shifted_distr <- function(distr, pzero){
  distr <- dist_transformed(distr, function(x) x + 1, function(x) x - 1)
  dist_inflated(distr, pzero, 0)
}

negative_mass_to_zero <- function(distr) {
  do.call(c, lapply(as.list(distr), \(d) {
    p0 <- cdf(d, 0)[[1]]
    dist_truncated(d, lower = 0) |> dist_inflated(prob = p0)
  }))
}

fit_nbinom <- function(y) {
  if (length(y) == 0 || all(y == 0)) {
    return(c(size = 100, prob = 1 - .STATICDISTR_EPSILON))
  }

  fit <- tryCatch(
    nloptr(
      x0 = c(max(mean(y), .STATICDISTR_EPSILON), 0.5),
      eval_f = function(x) -mean(dnbinom(y, x[1], x[2], log = TRUE)),
      lb = c(.STATICDISTR_EPSILON, .STATICDISTR_EPSILON),
      ub = c(Inf, 1 - .STATICDISTR_EPSILON),
      opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500)
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || is.null(fit$solution)) {
    mu <- mean(y)
    sigmasq <- var(y)
    if (!is.na(sigmasq) && sigmasq > mu + .STATICDISTR_EPSILON) {
      size <- (mu^2) / (sigmasq - mu)
    } else {
      size <- 100
    }
    prob <- min(size / (size + mu), 1 - .STATICDISTR_EPSILON)
    return(c(size = size, prob = prob))
  }

  c(size = fit$solution[1], prob = fit$solution[2])
}
