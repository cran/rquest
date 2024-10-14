#' qineq
#' @description
#' carry out hypothesis tests and obtain associated confidence intervals for quantile based inequality measures
#' @details
#'
#' This function performs hypothesis testing and calculates the corresponding confidence intervals for inequality measures based on quantiles.
#' The available options for quantile based measures in argument `measure` are shown below.
#'
#' * `QRI`: Quantile Ratio Index (Prendergast & Staudte, 2018). This is the default choice.
#' * `G2`:  Quantile variant of the Gini index (Prendergast & Staudte, 2016a).
#'
#' The default `var.method="qor"` is to estimate the probability density function directly using the lognormal Quantile Optimality Ratio (QOR)
#' for choosing a suitable bandwidth (Prendergast & Staudte,2016b). Alternatively, the variances can be
#' estimated by inverting a density estimator evaluated at the quantiles and this can be done using `var.method = "density"`.
#'
#' For more information and further examples, see Prendergast, Dedduwakumara & Staudte (2024)
#'
#' @param x a numeric vector of data values.
#' @param y an optional second vector of data values for two-sample testing.
#' @param J number of grid points
#' @param measure quantile based inequality measure to be used.  Either "QRI" (default) or "G2".(See details).
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param quantile.type argument for the quantile function.  Default is set to 8 so that output is consistent with default quantile function use and other functions such as IQR (see help file for `quantile()`
#' for more details)
#' @param var.method approach use to estimate the quantile density function.  Either "qor"(default) or "density".(See details).
#' @param conf.level coverage for the estimated confidence interval.
#' @param true.ineq the specified hypothesized value of the inequality measure or the difference of the inequality depending on whether it was a one-sample test or a two-sample test.
#' @return hypothesis test results and associated confidence interval (a list with class "htest")
#' @references
#'
#' Prendergast, L.A., & Staudte, R.G. (2016a). Quantile versions of the Lorenz curve. Electronic Journal of
#' Statistics, 10(2), 1896 – 1926.
#'
#' Prendergast, L. A., & Staudte, R. G. (2016b). Exploiting the quantile optimality ratio in finding confidence intervals for quantiles. Stat, 5(1), 70-81
#'
#' Prendergast, L. A., & Staudte, R. G. (2018). A simple and effective inequality measure. The American Statistician, 72(4), 328-343.
#'
#' Prendergast, L. A., Dedduwakumara, D.S. & Staudte, R.G. (2024) rquest: An R package for hypothesis tests and confidence intervals
#' for quantiles and summary measures based on quantiles, preprint, pages 1-13
#'
#' @export
#'
#' @examples
#' # Create some data
#' x <- c(8.43,7.08,8.79,8.88,7.87,5.94,8.79,5.46,8.11,7.08)
#' y <- c(13.44,13.65,14.77,9.51,14.07,10.92,11.59,13.42,8.93,10.88)
#'
#' # Two sample hypothesis test for the QRI measure
#' qineq(x,y)

qineq <- function(x, y = NULL, J = 100, measure = "QRI",
                  alternative = c("two.sided", "less", "greater"),
                  quantile.type = 8, var.method = "qor",
                  conf.level = 0.95, true.ineq = 0.5){

  alternative <- match.arg(alternative)
  if(is.null(y)){
    samples <- "One sample"
    dname <- deparse(substitute(x))
  } else{
    samples <- "Two sample"
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  }

  if(any(is.na(x))){
    count.x.na <- sum(is.na(x))
    warning(paste0(count.x.na), " missing values removed in ", deparse(substitute(x)), ".\n")
    x <- na.omit(x)
  }

  if(any(is.na(y))){
    count.y.na <- sum(is.na(y))
    warning(paste0(count.y.na), " missing values removed in ", deparse(substitute(y)), ".\n")
    y <- na.omit(y)
  }

  alpha <- 1 - conf.level
  crit <- qnorm(1 - alpha/2)

  p <- (1:J - 0.5)/J
  covQ <- qcov(x, sort(c(p/2, 1- p/2)), method = "qor", quantile.type = quantile.type)$cov # would be good to check
  if(measure == "QRI" | measure == "G2"){
    Rpx <- R(x, p, quantile.type = quantile.type)
    if(measure == "QRI"){
      multp <- rep(1, J)
      mult <- 1
    } else{
      multp <- p
      mult <- 2
    }
    estx <- (mult/J)*sum(multp - multp*Rpx)

    names(estx) <- measure
    method <- paste(samples, "test of the", measure)

    colR <- matrix(Rpx, J, J, byrow = FALSE)
    rowR <- matrix(Rpx, J, J, byrow = TRUE)
    covQtl <- covQ[1:J, 1:J]; covQtr <- covQ[1:J, (2*J):(J+1)]; covQbl <- covQ[(2*J):(J+1), 1:J]; covQbr <- covQ[(2*J):(J+1), (2*J):(J+1)]
    covR <- (multp %*% t(multp)) * ((1/quantile(x, 1 - p/2)) %*% t(1/quantile(x, 1 - p/2))) * (covQtl - rowR*covQtr - colR*covQbl + colR*rowR*covQbr)

    sterrx <- sqrt(mult^2*sum(covR)/J^2)

    if(!is.null(y)){
      Rpy <- R(y, p, quantile.type = quantile.type)
      esty <- (mult/J)*sum(multp - multp*Rpy)

      names(esty) <- measure

      colR <- matrix(Rpy, J, J, byrow = FALSE)
      rowR <- matrix(Rpy, J, J, byrow = TRUE)
      covQtl <- covQ[1:J, 1:J]; covQtr <- covQ[1:J, (2*J):(J+1)]; covQbl <- covQ[(2*J):(J+1), 1:J]; covQbr <- covQ[(2*J):(J+1), (2*J):(J+1)]
      covR <- (multp %*% t(multp)) * ((1/quantile(x, 1 - p/2)) %*% t(1/quantile(x, 1 - p/2))) * (covQtl - rowR*covQtr - colR*covQbl + colR*rowR*covQbr)

      sterry <- sqrt(mult^2*sum(covR)/J^2)
    }

  } else{
    stop("Only QRI and the G2 measures are supported at the moment.\n")
  }

  if (samples == "Two sample"){
    est <- estx - esty
    sterr <- sqrt(sterrx^2 + sterry^2)
    names(est) <- paste0("difference in ", measure)
    if(true.ineq == 0.5) true.ineq <- 0
  } else{
    est <- estx
    sterr <- sterrx
  }

  test.stat <- (est - true.ineq)/sterr
  names(test.stat) <- "Z"

  if(alternative == "less"){
    pval <- pnorm(test.stat)
    ci <- c(0, est + qnorm(conf.level)*sterr)
  } else if (alternative == "greater"){
    pval <- pnorm(test.stat, lower.tail = FALSE)
    ci <- c(est - qnorm(conf.level)*sterr, 1)
  } else{
    pval <- 2*(1 - pnorm(abs(test.stat)))
    ci <- est + c(-1, 1)*crit*sterr
  }

  attr(ci, 'conf.level') <- conf.level
  names(true.ineq) <- names(est)


  ineqres <- list(method = method,
                  data.name = dname,
                  statistic = test.stat,
                  parameter = NULL,
                  p.value = pval,
                  alternative = alternative,
                  estimate = est,
                  null.value = true.ineq,
                  conf.int = ci)

  class(ineqres) <- "htest"
  return(ineqres)
}
