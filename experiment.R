library(fable.intermittent)
library(fabletools)
library(tsibble)
library(dplyr)

### EXPERIMENT WITH AUTO DATASET ###

data("auto")
h <- 6

test_start <- auto$index |> 
  unique() |> 
  sort() |> 
  tail(h) |>
  head(1)

fit <- auto |>
  filter(index <= yearmonth("2011 Jun")) |>
  model(
    empdistr = EMPDISTR(value),
    staticdistr = STATICDISTR(value),
    betanbb = BETANBB(value),
    gampois = GAMPOISB(value),
    wss = WSS(value),
    vz = VZ(value),
    hspes = HSPES(value),
    negbines = NEGBINES(value),
    twees = TWEES(value),
    marwal = MARWAL(value)
  )

fc <- fit |>
  forecast(h = "6 months")

res <- fc |>
  accuracy(auto, measures = list(
    RMSSE = RMSSE, 
    pinball_loss = pinball_loss
    )) |>
  group_by(.model) |>
  summarise(RMSSE = mean(RMSSE), pinball_loss  = mean(pinball_loss))

print(res)

### EXPERIMENT WITH RAF DATASET ###

data("raf")
h <- 12

test_start <- raf$index |> 
  unique() |> 
  sort() |> 
  tail(h) |>
  head(1)

fit <- raf |>
  filter(index < test_start) |>
  model(
    empdistr = EMPDISTR(value),
    staticdistr = STATICDISTR(value),
    betanbb = BETANBB(value),
    gampois = GAMPOISB(value),
    wss = WSS(value),
    vz = VZ(value),
    hspes = HSPES(value),
    negbines = NEGBINES(value),
    twees = TWEES(value),
    marwal = MARWAL(value)
  )

fc <- fit |>
  forecast(h = h)

res <- fc |>
  accuracy(raf, measures = list(
    RMSSE = RMSSE, 
    pinball_loss = pinball_loss
  )) |>
  group_by(.model) |>
  summarise(RMSSE = mean(RMSSE), pinball_loss  = mean(pinball_loss))

print(res)
  