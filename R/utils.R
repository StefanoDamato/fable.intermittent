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