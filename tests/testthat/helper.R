# Test Data Helper for fable.intermittent
# Single tsibble with 5 different time series stacked, keyed by "series" variable
# Same test time series from probintermittent
library(tibble)
library(tsibble)

set.seed(42)
startdate <- as.Date("1997-08-12")

y1 <- c(rpois(30, 1.3), rpois(20, 0.3))
ts1 <- tsibble::tsibble(
  time = tsibble::yearmonth(startdate) + seq_len(length(y1)) - 1,
  value = y1,
  series = "TS1",
  index = "time",
  key = "series"
)

y2 <- rep(0, 20)
y2[sample.int(20, 1)] <- 1 + rpois(1, rgamma(1, 10, 1))
ts2 <- tsibble::tsibble(
  time = tsibble::yearquarter(startdate) + seq_len(length(y2)) - 1,
  value = y2,
  series = "TS2",
  index = "time",
  key = "series"
)

y3 <- rep(0, 150)
y3[sample.int(100, round(100 * runif(1, 0.1, 0.9)))] <- 1
y3[y3 == 1] <- rnbinom(sum(y3 == 1), 2, runif(1))
ts3 <- tsibble::tsibble(
  time = tsibble::yearweek(startdate) + seq_len(length(y3)) - 1,
  value = y3,
  series = "TS3",
  index = "time",
  key = "series"
)

y4 <- c(1, rep(0, 39))
ts4 <- tsibble::tsibble(
  time = startdate + seq_len(length(y4)) - 1,
  value = y4,
  series = "TS4",
  index = "time",
  key = "series"
)

y5 <- c(rpois(200, 2.5), rpois(200, 0.8), rpois(200, 1.5))
ts5 <- tsibble::tsibble(
  time = tsibble::yearweek(startdate) + seq_len(length(y5)) - 1,
  value = y5,
  series = "TS5",
  index = "time",
  key = "series"
)

y6 <- rnbinom(2000, 3 + 2.5 * sin(seq_len(2000) * 2* pi / 7), 0.6)
ts6 <- tsibble::tsibble(
  time = startdate + seq_len(length(y6)) - 1,
  value = y6,
  series = "TS6",
  index = "time",
  key = "series"
)

test_data = list(ts1, ts2, ts3, ts4, ts5, ts6) 
