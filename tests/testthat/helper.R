# Test Data Helper for fable.intermittent
# Single tsibble with 5 different time series stacked, keyed by "series" variable
# Same test time series from probintermittent
library(tibble)
library(tsibble)

set.seed(42)


ts1 <- tibble::tibble(
  time = 1:50,
  value = c(rpois(30, 1.3), rpois(20, 0.3)),
  series = "TS1"
)

y2 <- rep(0, 20)
y2[sample.int(20, 1)] <- 1 + rpois(1, rgamma(1, 10, 1))
ts2 <- tibble::tibble(
  time = 1:20,
  value = y2,
  series = "TS2"
)

y3 <- rep(0, 100)
y3[sample.int(100, round(100 * runif(1, 0.1, 0.9)))] <- 1
y3[y3 == 1] <- rnbinom(sum(y3 == 1), 2, runif(1))
ts3 <- tibble::tibble(
  time = 1:100,
  value = y3,
  series = "TS3"
)


ts4 <- tibble::tibble(
  time = 1:40,
  value = c(1, rep(0, 39)),
  series = "TS4"
)

ts5 <- tibble::tibble(
  time = 1:600,
  value = c(
    rpois(200, 2.5),  # High demand period
    rpois(200, 0.8),  # Medium demand period
    rpois(200, 1.5)   # Another medium period
  ),
  series = "TS5"
)

# Combine all into a single tsibble with key
test_data <- dplyr::bind_rows(ts1, ts2, ts3, ts4, ts5) |>
  tsibble::as_tsibble(index = time, key = series)
