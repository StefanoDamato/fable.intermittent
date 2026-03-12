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
