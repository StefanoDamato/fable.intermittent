#' RAF Spare Parts Demand Dataset
#'
#' A dataset of monthly demand for spare parts from the Royal Air Force (RAF).
#' The data contains 5000 intermittent time series, each spanning 84 monthly
#' periods from January 1996 to December 2002. This is a widely used benchmark
#' dataset for intermittent demand forecasting.
#'
#' @format A tsibble with 420,000 rows and 3 variables:
#' \describe{
#'   \item{series_id}{Character. Unique identifier for each time series.}
#'   \item{index}{Date (yearmonth). The monthly time index.}
#'   \item{value}{Numeric. The demand quantity for the given month.}
#' }
#'
#' @source Syntetos, A. A., & Boylan, J. E. (2005). The accuracy of
#'   intermittent demand estimates. \emph{International Journal of Forecasting},
#'   21(2), 303--314.
#'
#'   Available at
#'   \url{https://github.com/canerturkmen/gluon-ts/tree/intermittent-datasets/datasets}.
#'
#' @examples
#' library(tsibble)
#' raf
"raf"

#' Automotive Spare Parts Demand Dataset
#'
#' A dataset of monthly demand for automotive spare parts. The data contains
#' 3000 intermittent time series, each spanning 24 monthly periods from
#' January 2010 to December 2011.
#'
#' @format A tsibble with 72,000 rows and 3 variables:
#' \describe{
#'   \item{series_id}{Character. Unique identifier for each time series.}
#'   \item{index}{Date (yearmonth). The monthly time index.}
#'   \item{value}{Numeric. The demand quantity for the given month.}
#' }
#'
#' @source Turkmen, A. C., Januschowski, T., Wang, Y., & Cemgil, A. T. (2021).
#'   Forecasting intermittent and sparse time series: A unified probabilistic
#'   framework via deep renewal processes. \emph{PLOS ONE}, 16(11), e0259764.
#'
#'   Available at
#'   \url{https://github.com/canerturkmen/gluon-ts/tree/intermittent-datasets/datasets}.
#'
#' @examples
#' library(tsibble)
#' auto
"auto"

#' Pasta Sales Dataset
#'
#' Daily sales and promotional data for pasta products from an Italian grocery
#' store, spanning five years from January 2014 to December 2018. The dataset
#' covers 118 Stock Keeping Units (SKUs) organised across four brands (B1--B4),
#' making it a standard benchmark for hierarchical and intermittent demand
#' forecasting. Each SKU-day combination records the quantity sold and a binary
#' promotional indicator.
#'
#' @format A tsibble with 212,164 rows and 5 variables:
#' \describe{
#'   \item{week}{Date. Daily time index.}
#'   \item{brand}{Character. Brand identifier (B1, B2, B3, B4), with 42, 45,
#'     21, and 10 SKUs respectively.}
#'   \item{product}{Character. SKU number within the brand.}
#'   \item{value}{Numeric. Daily quantity sold.}
#'   \item{promotion}{Integer. Binary promotional indicator (1 = on promotion,
#'     0 = not on promotion).}
#' }
#'
#' @source Mancuso, P., Piccialli, V., & Sudoso, A. M. (2021).
#'   A machine learning approach for forecasting hierarchical time series.
#'   \emph{Expert Systems with Applications}, 182, 115102.
#'   \doi{10.1016/j.eswa.2021.115102}
#'
#'   Dataset available at
#'   \url{https://data.mendeley.com/datasets/njdkntcpc9/1}
#'   under a Creative Commons Attribution 4.0 International licence
#'   (CC BY 4.0).
#'
#' @examples
#' library(tsibble)
#' pasta
"pasta"
