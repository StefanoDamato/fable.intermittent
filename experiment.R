library(fable.intermittent)
library(fabletools)
library(tsibble)
library(dplyr)

### EXPERIMENT WITH AUTO DATASET ###

data("auto")

fit <- auto |>
  filter(index <= yearmonth("2011 Jun")) |>
  model(
    empdistr = EMPDISTR(value),
    staticdistr = STATICDISTR(value)
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

fit <- raf |>
  filter(index <= yearmonth("2001 Dec")) |>
  model(
    empdistr = EMPDISTR(value),
    staticdistr = STATICDISTR(value)
  )

fc <- fit |>
  forecast(h = "1 year")

res <- fc |>
  accuracy(raf, measures = list(
    RMSSE = RMSSE, 
    pinball_loss = pinball_loss
  )) |>
  group_by(.model) |>
  summarise(RMSSE = mean(RMSSE), pinball_loss  = mean(pinball_loss))

print(res)
  