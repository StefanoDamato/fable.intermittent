# Build package datasets from row-wise CSV files.
# Run from package root: source("data-raw/build_datasets.R")

build_index <- function(start, interval, n_periods) {
  offset <- 0:(n_periods - 1)

  switch(
    interval,
    month = tsibble::yearmonth(start) + offset,
    quarter = tsibble::yearquarter(start) + offset,
    week = as.Date(start) + (7 * offset),
    day = as.Date(start) + offset,
    year = seq.Date(as.Date(start), by = "year", length.out = n_periods),
    stop("Unsupported interval: ", interval)
  )
}

csv_to_tsibble <- function(csv_path, start, interval = "month", key_prefix = "TS") {
  interval <- match.arg(interval, c("month", "quarter", "week", "day", "year"))
  raw <- utils::read.csv(csv_path, check.names = FALSE)
  original_cols <- names(raw)
  n_periods <- length(original_cols)
  index_values <- build_index(start = start, interval = interval, n_periods = n_periods)

  out <- raw |>
    dplyr::mutate(series_id = paste0(key_prefix, dplyr::row_number())) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(original_cols),
      names_to = "period_col",
      values_to = "value"
    ) |>
    dplyr::mutate(
      period_index = match(period_col, original_cols),
      index = index_values[period_index],
      value = as.numeric(value)
    ) |>
    dplyr::select(series_id, index, value)

  tsibble::as_tsibble(out, key = series_id, index = index)
}

dataset_specs <- list(
  list(
    name = "auto",
    csv = file.path("data-raw", "Auto.csv"),
    start = "2010-01-01",
    interval = "month"
  ),
  list(
    name = "raf",
    csv = file.path("data-raw", "RAF.csv"),
    start = "1996-01-01",
    interval = "month"
  )
)

csv_paths <- vapply(dataset_specs, function(spec) spec$csv, character(1))
stopifnot(all(file.exists(csv_paths)))

datasets <- lapply(dataset_specs, function(spec) {
  csv_to_tsibble(
    csv_path = spec$csv,
    start = spec$start,
    interval = spec$interval,
  )
})
names(datasets) <- vapply(dataset_specs, function(spec) spec$name, character(1))

if (!dir.exists("data")) {
  dir.create("data", recursive = TRUE)
}

data_env <- list2env(datasets, parent = emptyenv())
for (obj_name in names(datasets)) {
  save(
    list = obj_name,
    file = file.path("data", paste0(obj_name, ".rda")),
    envir = data_env,
    compress = "xz"
  )
}

