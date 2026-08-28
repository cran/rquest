#' Hypothesis Tests and Confidence Intervals for the Robust Coefficient of Variation
#' @description
#' carry out hypothesis tests and obtain associated confidence intervals for robust versions of the coefficient of variation
#' @details
#'
#' This function performs hypothesis tests and constructs confidence intervals for robust versions of the coefficient of variation (rCV). It can be used for either one-sample or two-sample inference.
#'
#' The type of robust coefficient of variation is determined by the argument \code{numerator}. Two choices are available (Arachchige & Prendergast, 2022):
#'
#' \itemize{
#' \item \code{"mad"}: based on the median absolute deviation (MAD). This is the default and most widely used choice.
#' \item \code{"iqr"}: based on the interquartile range (IQR).
#' }
#'
#' Let \eqn{m} denote the sample median. Then the MAD-based robust coefficient of variation is
#' \deqn{\widehat{\mathrm{rCV}} = 1.4826 \times \frac{\mathrm{MAD}}{m}.}
#' The multiplier 1.4826 makes the estimator comparable to the usual coefficient of variation under a Gaussian model.
#'
#' If the IQR is preferred, then
#' \deqn{\widehat{\mathrm{rCV}} = 0.75 \times \frac{\mathrm{IQR}}{m}.}
#' The multiplier 0.75 is the corresponding scaling constant for the IQR-based version.
#'
#' For \code{numerator = "mad"}, Wald-type confidence intervals are computed using an asymptotic variance estimate for the rCV. The default is to perform the calculations on the log scale (\code{log.transf = TRUE}) and then, when \code{back.transf = TRUE}, exponentiate the estimate and interval endpoints back to the original rCV scale.
#'
#' For \code{numerator = "iqr"}, the rCV is treated as a ratio of linear combinations of quantiles and the function \code{q.test} is used to carry out the hypothesis test and compute the associated interval estimate. In this case, \code{var.method} controls how the quantile density is estimated.
#'
#' Additional arguments in \code{...} are passed to \code{density} when \code{numerator = "mad"} and to \code{q.test} when \code{numerator = "iqr"}.
#' @param x a numeric vector of data values.
#' @param y an optional second numeric vector of data values for two-sample testing.
#' @param numerator a character string specifying the numerator used in the robust coefficient of variation. Must be either \code{"mad"} or \code{"iqr"}.
#' @param alternative a character string specifying the alternative hypothesis; must be one of \code{"two.sided"} (default), \code{"less"}, or \code{"greater"}.
#' @param quantile.type argument for the quantile function. Default is \code{8} so that output is consistent with the default quantile function use and with other functions such as \code{IQR} (see \code{?quantile} for details).
#' @param var.method method used to estimate the quantile density function when \code{numerator = "iqr"}. Either \code{"qor"} (default) or \code{"density"}.
#' @param conf.level coverage level for the confidence interval.
#' @param true.rcv the hypothesized value of the robust coefficient of variation, or of the difference/ratio in the two-sample case, depending on the transformation used.
#' @param log.transf logical; if \code{TRUE}, the test and interval are computed on the log scale.
#' @param back.transf logical; if \code{TRUE}, estimates and confidence intervals are back-transformed to the original rCV scale using \code{exp}.
#' @param ... additional arguments to be passed to \code{density} when \code{numerator = "mad"} and to \code{q.test} when \code{numerator = "iqr"}.
#' @return an object of class \code{"htest"} containing the test result, estimate, null value, p-value, and confidence interval.
#' @references
#' Arachchige, C. N. P. G., & Prendergast, L. A. (2022). Robust analogs to the coefficient of variation. \emph{Journal of Applied Statistics}, \strong{49}(2), 268--290.
#' @seealso
#' \code{\link{qcov}} for covariance estimation of sample quantiles,
#' \code{\link{qrcov}} for covariance estimation of ratios of linear combinations of quantiles, and
#' \code{\link{q.test}} for hypothesis tests and confidence intervals based on ratios of linear combinations of quantiles.
#' @export
#' @examples
#' set.seed(123)
#' x1 <- rnorm(100, 8)
#' x2 <- rnorm(120, 10)
#'
#' # One-sample tests
#' rcv.test(x1)
#' rcv.test(x1, numerator = "iqr")
#'
#' # Two-sample tests
#' rcv.test(x1, x2)
#' rcv.test(x1, x2, numerator = "iqr")

rcv.test <- function (x, y = NULL, numerator = "mad", quantile.type = 8, var.method = "qor",
          alternative = c("two.sided", "less", "greater"), conf.level = 0.95,
          true.rcv = 0, log.transf = TRUE, back.transf = TRUE, ...)
{
  if (!is.numeric(x))
    stop("Argument 'x' must be numeric.")
  if (is.null(y)) {
    samples <- "One sample"
    dname <- deparse(substitute(x))
  }
  else {
    if (!is.numeric(y))
      stop("argument 'y' must be numeric")
    samples <- "Two sample"
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  }
  if (anyNA(x)) {
    count.x.na <- sum(is.na(x))
    warning(paste0(count.x.na), " missing values removed in ",
            deparse(substitute(x)), ".\n")
    x <- na.omit(x)
  }
  if (anyNA(y)) {
    count.y.na <- sum(is.na(y))
    warning(paste0(count.y.na), " missing values removed in ",
            deparse(substitute(y)), ".\n")
    y <- na.omit(y)
  }
  alternative <- match.arg(alternative)

  nx <- length(x)
  ny <- length(y)
  alpha <- 1 - conf.level
  crit <- qnorm(1 - alpha/2)

  if (numerator == "iqr"){
    out <- q.test(x = x, y = y, measure = "rCViqr", quantile.type = quantile.type,
                  true.q = true.rcv, log.transf = log.transf, back.transf = back.transf,
                  conf.level = conf.level, var.method = var.method, alternative = alternative, ...)
  } else if(numerator == "mad"){
    densx <- density(x, ...)
    med.estx <- median(x)
    mad.estx <- mad(x, constant = 1)
    estx <- 1.4826*mad.estx/med.estx

    fFinvx <- approx(densx$x, densx$y, xout = med.estx + c(-1, 1, 0)*mad.estx)$y
    F.estx <- ecdf(x)
    FFinvx <- F.estx(med.estx + c(-1, 1)*mad.estx)

    C1x <- fFinvx[1] + fFinvx[2]
    C3x <- fFinvx[1] - fFinvx[2]
    C2x <- C3x^2 + 4*C3x*fFinvx[3]*(1 - FFinvx[2] - FFinvx[1])

    rho1x <- 1/fFinvx[3]^2/4
    rho2x <- (1 + C2x/fFinvx[3]^2)/4/C1x^2
    rho12x <- (1 - 4*FFinvx[1] + C3x/fFinvx[3])/4/C1x^2

    sterrx <- sqrt(estx^2*(rho1x/med.estx^2 +rho2x/mad.estx^2 - 2*rho12x/med.estx/mad.estx)/nx)

    if (log.transf) {
      if (estx <= 0) {
        stop("Estimates must be positive to use the log transformation.\n")
      }
      sterrx <- sterrx/estx
      estx <- log(estx)
      if (back.transf)
        transf.text <- NULL
      else transf.text <- "(log transformed)"
    }
    else transf.text <- NULL

    if (!is.null(y)) {
      densy <- density(y, ...)
      med.esty <- median(y)
      mad.esty <- mad(y, constant = 1)
      esty <- 1.4826*mad.esty/med.esty

      fFinvy <- approx(densy$x, densy$y, xout = med.esty + c(-1, 1, 0)*mad.esty)$y
      F.esty <- ecdf(y)
      FFinvy <- F.esty(med.esty + c(-1, 1)*mad.esty)

      C1y <- fFinvy[1] + fFinvy[2]
      C3y <- fFinvy[1] - fFinvy[2]
      C2y <- C3y^2 + 4*C3y*fFinvy[3]*(1 - FFinvy[2] - FFinvy[1])

      rho1y <- 1/fFinvy[3]^2/4
      rho2y <- (1 + C2y/fFinvy[3]^2)/4/C1y^2
      rho12y <- (1 - 4*FFinvy[1] + C3y/fFinvy[3])/4/C1y^2

      sterry <- sqrt(esty^2*(rho1y/med.esty^2 +rho2y/mad.esty^2 - 2*rho12y/med.esty/mad.esty)/ny)

      if (log.transf) {
        if (esty <= 0) {
          stop("Estimates must be positive to use the log transformation.\n")
        }
        sterry <- sterry/esty
        esty <- log(esty)
      }
    }
    measure.name <- "Robust CV"
    method <- paste(samples, "test of the robust coefficient of variation (MAD/median)")

    if (samples == "Two sample") {
      est <- estx - esty
      sterr <- sqrt(sterrx^2 + sterry^2)
    }
    else {
      est <- estx
      sterr <- sterrx
    }
    test.stat <- (est - true.rcv)/sterr
    names(test.stat) <- "Z"
    if (alternative == "less") {
      pval <- pnorm(test.stat)
      ci <- c(-Inf, est + qnorm(conf.level) * sterr)
    }
    else if (alternative == "greater") {
      pval <- pnorm(test.stat, lower.tail = FALSE)
      ci <- c(est - qnorm(conf.level) * sterr, Inf)
    }
    else {
      pval <- 2 * (1 - pnorm(abs(test.stat)))
      ci <- est + c(-1, 1) * crit * sterr
    }
    if (log.transf & back.transf) {
      est <- exp(est)
      ci <- exp(ci)
      if (true.rcv == 0)
        true.rcv <- 1
    }
    attr(ci, "conf.level") <- conf.level
    names(est) <- "Robust CV"
    if (samples == "Two sample") {
      names(est) <- paste0(names(est), "s")
      if (log.transf & back.transf)
        names(est) <- paste("ratio of", names(est))
      else names(est) <- paste("difference in", names(est))
    }
    names(est) <- paste(names(est), transf.text)
    names(true.rcv) <- names(est)
    out <- list(method = method, data.name = dname, statistic = test.stat,
                 parameter = NULL, p.value = pval, alternative = alternative,
                 estimate = est, null.value = true.rcv, conf.int = ci)
    class(out) <- "htest"

  } else stop("Argument 'numerator' must be either 'mad' or 'iqr'.\n")


  return(out)
}
