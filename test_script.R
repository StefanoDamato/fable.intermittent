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
ts <- test_data[[5]]
#ts1 <- raf |>
  #dplyr::filter(series_id == "TS1")
fit <- ts |> model(GAMPOISB(value))
fc <- fit |> forecast(h = 10)
fc$.mean
gen <- fit |> generate(h = 10)
# =======================================================================



# example ======================================================================
undebug(train_betanbb)
undebug(forecast.BETANBB)
h <- 10

model_builders <- list(
  gampoisb = function() GAMPOISB(value),
  betanbb = function() BETANBB(value),
  empdistr = function() EMPDISTR(value),
  staticdistr = function() STATICDISTR(value),
  negbines = function() NEGBINES(value),
  hspes = function() HSPES(value),
  wss = function() WSS(value),
  vz = function() VZ(value),
  marwal = function() MARWAL(value)
  #poisgas = function() POISGAS(value)
)

run_models_one_series <- function(ts_data, h, model_builders) {
  series_id <- unique(ts_data$series)
  if (length(series_id) != 1) {
    stop("Each element of test_data must contain exactly one series.")
  }

  fc_list <- list()
  fit_list <- list()

  for (model_name in names(model_builders)) {
    fit_obj <- tryCatch(
      ts_data |> model(.tmp = model_builders[[model_name]]()),
      error = function(e) {
        message("[", series_id, "] model '", model_name, "' failed in training: ", e$message)
        NULL
      }
    )

    if (is.null(fit_obj)) {
      next
    }

    fc_obj <- tryCatch(
      fit_obj |> forecast(h = h),
      error = function(e) {
        message("[", series_id, "] model '", model_name, "' failed in forecast: ", e$message)
        NULL
      }
    )

    fit_vals <- tryCatch(
      fit_obj |> fitted(),
      error = function(e) {
        message("[", series_id, "] model '", model_name, "' failed in fitted(): ", e$message)
        NULL
      }
    )

    if (!is.null(fc_obj)) {
      fc_list[[model_name]] <- dplyr::mutate(fc_obj, .model = model_name)
    }
    if (!is.null(fit_vals)) {
      fit_list[[model_name]] <- dplyr::mutate(fit_vals, .model = model_name)
    }
  }

  if (length(fc_list) == 0) {
    message("[", series_id, "] no model produced forecasts.")
    return(NULL)
  }

  fc_all <- dplyr::bind_rows(fc_list)
  fit_all <- dplyr::bind_rows(fit_list)

  p <- fc_all |>
    autoplot(ts_data, level = NULL) +
    geom_line(data = fit_all, aes(x = time, y = .fitted, colour = .model), alpha = 0.6) +
    labs(
      title = paste0("Forecast comparison - ", series_id),
      x = "Time",
      y = "Value",
      colour = "Model"
    )

  list(series = series_id, forecast = fc_all, fitted = fit_all, plot = p)
}

results <- lapply(test_data, run_models_one_series, h = h, model_builders = model_builders)
results <- Filter(Negate(is.null), results)

for (res in results) {
  print(res$plot)
}

# Combined plot across all series =============================================
if (length(results) > 0) {
  all_fc <- dplyr::bind_rows(lapply(results, function(r) r$forecast))
  all_fit <- dplyr::bind_rows(lapply(results, function(r) r$fitted))
  all_orig <- dplyr::bind_rows(test_data)
  
  p_combined <- all_fc |>
    autoplot(all_orig, level = NULL) +
    geom_line(data = all_fit, aes(x = time, y = .fitted, colour = .model), alpha = 0.6) +
    facet_wrap(~series, scales = "free_y", ncol = 2) +
    labs(
      title = "All Models - All Series Combined",
      x = "Time",
      y = "Value",
      colour = "Model"
    ) +
    theme(legend.position = "bottom")
  
  print(p_combined)
}


Rcpp::sourceCpp("src/gas.cpp")


y <- ts1$value
psi0 <- 0.1
phi <- -.2
rho <- 0.15
xi0 <- -0.1
k <- 0.2
period <- 7

r_filter <- gas_filter(y, psi0, phi, rho, xi0, k, period = 5, distr = "pois")
cpp_filter <- gasFilter(y, psi0, phi, rho, xi0, k, period = 7)

r_filter$f
cpp_filter$f


gasFilterPois(y, psi0, phi, rho, xi0, k, period = 4)

