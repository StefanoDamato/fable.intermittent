test_that("raf is a well-formed tsibble", {
  expect_s3_class(raf, "tbl_ts")
  expect_equal(tsibble::key_vars(raf), "series_id")
  expect_equal(tsibble::index_var(raf), "index")
  expect_equal(tsibble::n_keys(raf), 5000)
  expect_false(any(tsibble::has_gaps(raf)$.gaps))
})

test_that("raf has the documented shape and no missing values", {
  expect_equal(dim(raf), c(420000, 3))
  expect_equal(names(raf), c("series_id", "index", "value"))
  expect_equal(range(raf$index), tsibble::yearmonth(c("1996 Jan", "2002 Dec")))
  expect_false(anyNA(raf$value))
  expect_true(all(raf$value >= 0))
})

test_that("auto is a well-formed tsibble", {
  expect_s3_class(auto, "tbl_ts")
  expect_equal(tsibble::key_vars(auto), "series_id")
  expect_equal(tsibble::index_var(auto), "index")
  expect_equal(tsibble::n_keys(auto), 3000)
  expect_false(any(tsibble::has_gaps(auto)$.gaps))
})

test_that("auto has the documented shape and no missing values", {
  expect_equal(dim(auto), c(72000, 3))
  expect_equal(names(auto), c("series_id", "index", "value"))
  expect_equal(range(auto$index), tsibble::yearmonth(c("2010 Jan", "2011 Dec")))
  expect_false(anyNA(auto$value))
  expect_true(all(auto$value >= 0))
})

test_that("pasta is a well-formed tsibble", {
  expect_s3_class(pasta, "tbl_ts")
  expect_equal(tsibble::key_vars(pasta), c("brand", "product"))
  expect_equal(tsibble::index_var(pasta), "index")
  expect_equal(tsibble::n_keys(pasta), 118)
  expect_false(any(tsibble::has_gaps(pasta)$.gaps))
})

test_that("pasta has the documented shape, SKU counts, and no missing values", {
  expect_equal(dim(pasta), c(215350, 5))
  expect_equal(names(pasta), c("index", "brand", "product", "value", "promotion"))
  expect_equal(range(pasta$index), as.Date(c("2014-01-02", "2018-12-31")))

  sku_counts <- tsibble::as_tibble(pasta) |>
    dplyr::distinct(brand, product) |>
    dplyr::count(brand) |>
    dplyr::pull(n, name = brand)
  expect_equal(sku_counts[c("B1", "B2", "B3", "B4")], c(B1 = 42, B2 = 45, B3 = 21, B4 = 10))

  expect_false(anyNA(pasta$value))
  expect_false(anyNA(pasta$promotion))
  expect_true(all(pasta$value >= 0))
  expect_true(all(pasta$promotion %in% c(0, 1)))
})

test_that("tinyM5 is a well-formed tsibble", {
  expect_s3_class(tinyM5, "tbl_ts")
  expect_equal(tsibble::key_vars(tinyM5), c("item_id", "store_id"))
  expect_equal(tsibble::index_var(tinyM5), "date")
  expect_equal(tsibble::n_keys(tinyM5), 280)
  expect_false(any(tsibble::has_gaps(tinyM5)$.gaps))
})

test_that("tinyM5 has the documented shape and column types", {
  expect_equal(dim(tinyM5), c(535640, 18))
  expect_equal(
    names(tinyM5),
    c("item_id", "dept_id", "cat_id", "store_id", "state_id", "value",
      "date", "wm_yr_wk", "weekday", "wday", "month", "year",
      "event_name_1", "event_type_1", "event_name_2", "event_type_2",
      "snap", "sell_price")
  )
  expect_equal(range(tinyM5$date), as.Date(c("2011-01-29", "2016-04-24")))
})

test_that("tinyM5 key and demand columns have no missing values", {
  key_and_demand_cols <- c("item_id", "store_id", "date", "value")
  expect_all_true(!sapply(as.data.frame(tinyM5)[key_and_demand_cols], anyNA))
})

test_that("tinyM5 demand values are non-negative", {
  expect_true(all(tinyM5$value >= 0))
})
