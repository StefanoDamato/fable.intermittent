#' Generate Intermittent Time Series Data
#'
#' Creates synthetic intermittent time series data in tsibble format.
#' Intermittent series are characterized by many zero values with sporadic non-zero observations.
#'
#' @param num_ts Number of time series to generate
#' @param len_ts Length of each time series
#' @param freq Frequency of the time series: "W" (weekly), "D" (daily), or "M" (monthly)
#' @param seasonality_w Weight of seasonality component (0 to 1, where 0 = no seasonality)
#' @param prob_nonzero Probability of a non-zero observation (between 0 and 1)
#' @param demand_size Shape parameter for gamma distribution (higher = lower variability)
#'
#' @return A tsibble object with columns: series_id, index, value (integer)
#'
#' @examples
#' # Generate 3 monthly intermittent series of length 100
#' ts_data <- generate_intermittent_ts(num_ts = 3, len_ts = 100, freq = "M", seasonality_w = 0.3)
#'
#' @export
generate_intermittent_ts <- function(num_ts, len_ts, freq = "M", seasonality_w = 0.2, prob_nonzero = 0.3, demand_size = 3) {
  # Validate inputs
  stopifnot(
    is.numeric(num_ts) && num_ts > 0,
    is.numeric(len_ts) && len_ts > 0,
    freq %in% c("W", "D", "M"),
    is.numeric(seasonality_w) && seasonality_w >= 0 && seasonality_w <= 1
  )
  
  # Generate dates based on frequency
  start_date <- as.Date("2026-01-01")
  
  index_values <- switch(
    freq,
    "W" = start_date + (0:(len_ts - 1)) * 7,
    "D" = start_date + (0:(len_ts - 1)),
    "M" = tsibble::yearmonth(start_date) + (0:(len_ts - 1))
  )
  
  # Generate time series data
  data_list <- lapply(1:num_ts, function(ts_id) {
    
    # Seasonal pattern
    if (freq == "M") {
      seasonal_period <- 12
    } else if (freq == "W") {
      seasonal_period <- 52
    } else {  # daily
      seasonal_period <- 365
    }
    
    seasonal <- seasonality_w * sin(2 * pi * (1:len_ts) / seasonal_period)
    
    # Generate intermittent values from gamma distribution
    # Sparse observations with zero-inflation
    sparse_indicator <- rbinom(len_ts, size = 1, prob = prob_nonzero)
    
    # Sample from gamma distribution with specified demand_size (shape parameter)
    # Scale is set to 1 so mean ≈ demand_size
    gamma_samples <- rgamma(len_ts, shape = demand_size, scale = 1)
    
    # Combine with seasonality and convert to integers
    values <- sparse_indicator * as.integer(gamma_samples * (1 + seasonal))
    
    data.frame(
      series_id = paste0("TS", ts_id),
      index = index_values,
      value = values,
      stringsAsFactors = FALSE
    )
  })
  
  # Combine all series
  data_combined <- do.call(rbind, data_list)
  rownames(data_combined) <- NULL
  
  # Convert to tsibble
  if (freq == "M") {
    ts_obj <- tsibble::as_tsibble(
      data_combined,
      key = series_id,
      index = index
    )
  } else if (freq == "W") {
    data_combined$index <- tsibble::yearweek(data_combined$index)
    ts_obj <- tsibble::as_tsibble(
      data_combined,
      key = series_id,
      index = index
    )
  } else {  # daily
    ts_obj <- tsibble::as_tsibble(
      data_combined,
      key = series_id,
      index = index
    )
  }
  
  return(ts_obj)
}
