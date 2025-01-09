#' R
#'
#' calculates the ratio of symmetric quantiles.
#'
#' @param x a numeric vector of data values.
#' @param p a numeric vector of probability values in \[0, 1\] indicating all quantiles to be estimated.
#' @param quantile.type argument for the quantile function.  Default is set to 8 so that output is consistent with default quantile function use and other functions such as IQR.
#' @return the ratio of symmetric quantiles
#' @noRd

R <- function(x, p, quantile.type = 8){
  xp2 <- quantile(x, p/2, quantile.type = quantile.type)
  x1p2 <- quantile(x, 1 - p/2, quantile.type = quantile.type)
  xp2/x1p2
}
