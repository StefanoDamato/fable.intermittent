#!/usr/bin/env Rscript
# =============================================================================
# Test Script for fable.intermittent
# Tests both fitting and prediction for WSS, VZ, and NEGBIN models
# =============================================================================

library(fable.intermittent)
library(fabletools)
library(tsibble)
library(dplyr)
library(ggplot2)

# Load test data
source("tests/testthat/helper.R")
devtools::load_all()


# debugger =========================================================
debug(train_vz)
debug(forecast.VZ)
ts1 <- test_data |> 
  dplyr::filter(series == "TS1")
fit <- ts1 |> model(VZ(value))
fc <- fit |> forecast(h = 10)
gen <- fit |> generate(h = 10)
# =======================================================================



# example ======================================================================
fit <- test_data |> 
  model(
    gampoisb = GAMPOISB(value),
    betanbb = BETANBB(value),
    empdistr = EMPDISTR(value),
    #staticdistr = STATICDISTR(value),
    negbines = NEGBINES(value),
    hspes = HSPES(value),
    wss = WSS(value),
    vz = VZ(value),
    )
fc <- fit |> forecast(h = 10)
fc <- fc |> dplyr::mutate(q0.9 = quantile(value, 0.9))

fits <- fit |> fitted()


fc |> autoplot(test_data, level=NULL)+
  geom_line(data = fits, aes(x = time, y = .fitted, colour = .model), alpha = 0.6) +
  geom_line(data = fc, aes(x = time, y = q0.9, colour = .model), linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~series, ncol=1, scales = 'free') +
  labs(
    x = "Time",
    y = "Value"
  )

test_data |> dplyr::filter(series == "TS2")
