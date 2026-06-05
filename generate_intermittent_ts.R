#' Generate Intermittent Time Series Data
#'
#' Creates synthetic intermittent time series data in tsibble format.
#' The occurrence is generated using a Bernoulli distribution, whose probability
#' can be influenced by a seasonal pattern. The demand size is generated using 
#' a gamma distribution, which can also be influenced by seasonality.
#'
#' @param num_ts Number of time series to generate
#' @param len_ts Length of each time series
#' @param freq Frequency of the time series: "W" (weekly), "D" (daily), or "M" (monthly)
#' @param prob_occurrence Probability of a non-zero observation (between 0 and 1)
#' @param occurrence_seasonal_w Weight of seasonality in occurrence (0 to 1, where 0 = no seasonality)
#' @param demand_size Mean and variance of the gamma distribution (higher = lower variability)
#' @param demand_seasonality_w Weight of seasonality component (0 to 1, where 0 = no seasonality)
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
generate_intermittent_ts <- function(num_ts, len_ts, freq = "M", 
                                     prob_occurrence = 0.3, occurrence_seasonal_w = 0.2,
                                     demand_size = 3, demand_seasonal_w = 0.2) {
  
  # Validate inputs
  stopifnot(
    is.numeric(num_ts) && num_ts > 0,
    is.numeric(len_ts) && len_ts > 0,
    freq %in% c("W", "D", "M"),
    is.numeric(prob_occurrence) && prob_occurrence >= 0 && prob_occurrence <= 1,
    is.numeric(occurrence_seasonal_w) && occurrence_seasonal_w >= 0 && occurrence_seasonal_w <= 1,
    is.numeric(demand_size) && demand_size > 0,
    is.numeric(demand_seasonal_w) && demand_seasonal_w >= 0 && demand_seasonal_w <= 1
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
    } else if (freq == "D") {  # daily
      seasonal_period <- 7
    }
    
    # Generate the seasonal weight
    seasonal_occurrence <- occurrence_seasonal_w * sin(2 * pi * (1:len_ts) / seasonal_period)
    seasonal_demand <- demand_seasonal_w * sin(2 * pi * (1:len_ts) / seasonal_period)
    
    # Generate the occurrence with a Bernoulli distribution 
    prob <- prob_occurrence + seasonal_occurrence*(0.5 + abs(0.5 - prob_occurrence))
    occurrence <- rbinom(len_ts, size = 1, prob = prob)
    
    # Generate the seasonal demand size using a gamma distribution
    mean_size <- demand_size * (1 + seasonal_demand)
    demand_size <- rgamma(len_ts, mean_size, 1)
    
    # Combine with seasonality and convert to integers
    values <- occurrence * ceiling(demand_size)
    
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


dataset <- generate_intermittent_ts(num_ts = 1, len_ts = 100, freq = "M")



