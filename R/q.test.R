#' q.test
#' @description
#' carry out hypothesis tests and obtain associated confidence intervals for linear combinations of quantiles,
#' and ratios of such linear combinations.
#' @details
#'
#' This function `q.test` performs hypothesis tests and calculates confidence intervals for linear combinations of quantiles.
#' The quantile measures that can be estimated are specified in the `measure` argument and are listed below.
#'
#' * `median`: The default choice.
#' * `iqr`: The interquartile range.
#' * `rCViqr`: The robust Coefficient of Variation measure (Arachchige et al.,2022)
#' * `bowley`, `groenR` and `groenL`: These choices are for Bowley's skew coefficient,
#' for the generalized measure, and the right and left skew measures by Groeneveld and Meeden (Groeneveld & Meeden, 1984, 2009). For each of these measures, the user may specify the choice of p using argument `p`.  If this is not used, then the default is `p=0.25`.
#' * `moors`: This Moors kurtosis measure (Moors, 1988).
#' * `lqw` and `rqw`: These are for the robust left and right tail weights (Brys et al., 2006).
#' * `qrxxyy`: A character string consisting of the first two characters "qr" and followed by four numbers will request a ratio of dependents quantiles (i.e. different quantiles from the same sample).The first two number digits (in place of "xx") will indicate the quantile for the numerator, and the second two numerical digits "yy" for the denominator.  For example, `qr9010` will estimate and test the ratio Q(0.9)/Q(0.1).
#'
#' The default `var.method="qor"` is to estimate the probability density function directly using the lognormal Quantile Optimality Ratio (QOR)
#' for choosing a suitable bandwidth (Prendergast & Staudte,2016). Alternatively, the variances can be
#' estimated by inverting a density estimator evaluated at the quantiles and this can be done using `var.method = "density"`.
#'
#' Additional to using a text string with argument measure to indicate which quantile-based measure is to be used (of those included for the `q.test` function),
#' users can also define and request their own.  For example, while the median is the default for a single quantile, other quantiles can also be requested
#' (e.g., `u = 0.25` will request the first quartile).  This option is also included for linear combinations of quantiles and ratios of linear combinations,
#' which can be done in two ways.  For a single linear combination, argument `u` is the vector of probability values defining the quantiles to be used,
#' and argument `coef` is a vector of coefficients for the linear combination.  If a ratio is needed, then `u` and `coef` are used for the numerator and `u2` and `coef2`
#' for the denominator linear combination.  Alternatively, a single `u` can be used to identify all quantiles for the ratio, and then `coef` can be a matrix (2 rows)
#' whose first row specifies the coefficients for the numerator and the second row for the denominator. For more information and further examples,
#' see Prendergast, Dedduwakumara & Staudte (2024), and the example code below shows how to obtain results for
#' the robust CVs using all three approaches (to achieve identical results).
#'
#' Further details on confidence intervals, including coverage properties of, can be found for ratios of dependent quantiles (Prendergast & Staudte, 2017),
#' ratios of independent quantiles and IQRs (Arachchige et al., 2021) and for robust CVs (Arachchige et al., 2022).
#'
#' @param x a numeric vector of data values.
#' @param y an optional second vector of data values for two-sample testing.
#' @param measure a character string specifying the quantile measure to be estimated (See details).
#' @param u a numeric vector of probability values in \[0, 1\] indicating all quantiles to be estimated.
#' @param coef a vector or matrix with two rows specifying the coefficients that define the linear combinations (coefficients must match the corresponding probability values in u). If coef is a vector then a single linear combination (LC) is computed. If it is a matrix, then first row defines the numerator LS and the second the denominator LC.
#' @param u2 a numeric vector of probability values in \[0, 1\] indicating all quantiles to be estimated for the denominator.
#' @param coef2 a vector specifying the coefficients that define the linear combination for the denominator.  This is can be used as an alternative to defining coef as matrix for ratios of linear combinations.
#' @param quantile.type argument for the quantile function.  Default is set to 8 so that output is consistent with default quantile function use and other functions such as IQR (see help file for `quantile()`
#' for more details)
#' @param var.method approach use to estimate the quantile density function.  Either "qor"(default) or "density" (See details).
#' @param alternative a character string for alternative hypothesis equal to one of "two.sided", "greater" or "less".
#' @param conf.level coverage for the estimated confidence interval.
#' @param true.q a numeric value for the true value under the null hypothesis test.
#' @param log.transf boolean indicating whether the a log transformation of the measure is to be used (i.e., estimates of the log of the measure are computed).
#' @param back.transf boolean indicating whether the measure and estimates should be back-transformed  to the original scale using exp.
#' @param min.q the lower bound for a one-sided confidence interval when alternative  argument if "less".
#' @param p optional value in (0, 1) for Bowley's generalized skewness coefficient.
#' @return hypothesis test results and associated confidence interval (a list with class "htest")
#' @references
#'
#' Arachchige, C. N., Cairns, M., & Prendergast, L. A. (2021). Interval estimators for ratios of independent quantiles and interquantile ranges.
#' Communications in Statistics-Simulation and Computation, 50(12), 3914-3930.
#'
#' Arachchige, C.N.P.G., Prendergast, L.A., & Staudte, R.G. (2022). Robust Analogs to the Coefficient of
#' Variation. Journal of Applied Statistics, 49(2), 268–290.
#'
#' Brys, G., Hubert, M., & Struyf, A. (2006). Robust measures of tail weight. Computational Statistics & Data
#' Analysis, 50(3), 733–759.
#'
#' Groeneveld, R. A., & Meeden, G. (1984). Measuring skewness and kurtosis. Journal of the Royal Statistical
#' Society Series D: The Statistician, 33(4), 391–399.
#'
#' Groeneveld, R. A., & Meeden, G. (2009). An improved skewness measure. Metron, 67(3), 325.
#'
#' Hyndman, R.J., & Fan, Y. (1996). Sample quantiles in statistical packages. The American Statistician, 50(4),
#' 361–365.
#'
#' Moors, J. J. A. (1988). A quantile alternative for kurtosis. Journal of the Royal Statistical Society:
#' Series D (The Statistician), 37(1), 25–32.
#'
#' Prendergast, L. A., & Staudte, R. G. (2016). Exploiting the quantile optimality ratio in finding confidence intervals for quantiles. Stat, 5(1), 70-81
#'
#' Prendergast, L. A., & Staudte, R. G. (2017). When large n is not enough–distribution-free interval estimators for ratios of quantiles.
#' The Journal of Economic Inequality, 15, 277-293.
#'
#' Prendergast, L. A., Dedduwakumara, D.S. & Staudte, R.G. (2024) rquest: An R package for hypothesis tests and confidence intervals
#' for quantiles and summary measures based on quantiles, preprint, pages 1-13
#'
#' @export
#'
#' @examples
#'
#' #  Create some data
#' x <- c(8.43,7.08,8.79,8.88,7.87,5.94,8.79,5.46,8.11,7.08)
#' y <- c(13.44,13.65,14.77,9.51,14.07,10.92,11.59,13.42,8.93,10.88)
#'
#' # One sample hypothesis test for the IQR
#' q.test(x, measure = "iqr")
#'
#' # Two samples hypothesis test for robust coefficient variations (0.75*IQR/median)
#' # with log transformation and back-transformation to the ratio scale,.
#' q.test(x, y, measure = "rCViqr", log.transf = TRUE, back.transf = TRUE)
#'
#' # The same two samples hypothesis test for robust coefficient variations (0.75*IQR/median)
#' # by using 'u',''u2','coef' and 'coef2' arguments.
#' u<-c(0.25,0.75)
#' coef<-0.75*c(-1,1)
#' u2<-0.5
#' coef2<-1
#' q.test(x,y,u=u,u2=u2,coef=coef,coef2=coef2,log.transf=TRUE,back.transf=TRUE)
#'
#' # The same two samples hypothesis test for robust coefficient variations (0.75*IQR/median)
#' # by using only 'u' and 'coef' arguments.
#' u<-c(0.25,0.5,0.75)
#' num <- 0.75*c(-1,0,1)
#' den <- c(0,1,0)
#' coef <- rbind(num, den)
#' q.test(x,y,u=u,coef=coef,log.transf=TRUE,back.transf=TRUE)


# Linear combination of quantiles -----------------------------------------

q.test <- function(x, y = NULL, measure = "median", u = NULL, coef = NULL, u2 = NULL, coef2 = NULL, quantile.type = 8, var.method = "qor",
                   alternative = c("two.sided", "less", "greater"),
                   conf.level = 0.95, true.q = 0, log.transf = FALSE, back.transf = FALSE, min.q = -Inf,
                   p = NULL){

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

  alternative <- match.arg(alternative)


  if(is.null(coef) & !is.null(u)) coef <- rep(1, length(u))

  if(!is.null(u2) | !is.null(coef2)){
    if(is.null(coef2)) coef2 <- rep(1, length(u2))
    if(is.null(u) | is.null(coef)){
      stop("When using u2 and coef2, you also need to specify both u and coef.\n")
    }
    if(is.null(u2)) stop("When using coef2 you also need to specify u2.\n")
    l <- length(coef)
    l2 <- length(u2)
    u <- c(u, u2)
    coef <- c(coef, rep(0, l2))
    coef2 <- c(rep(0, l), coef2)
    coef <- rbind(coef, coef2)
  }

  if(is.null(u) & is.null(coef)){
    if (measure == "rCViqr") {
      u <- c(0.25, 0.5, 0.75)
      coef <- matrix(c(-3/4, 0, 3/4, 0, 1, 0), nrow = 2, ncol = 3, byrow = TRUE)
      measure.name <- "Robust CV"
      method <- paste(samples, "test of the robust coefficient of variation (0.75*IQR/median)")
    } else if (measure == "iqr" | measure == "IQR") {
      u <- c(0.25, 0.75)
      coef <- c(-1, 1)
      measure.name <- "IQR"
      method <- paste(samples, "test of the interquartile range (IQR)")
    } else if (measure == "median" | measure == "med") {
      u <- 0.5
      coef <- 1
      measure.name <- "median"
      method <- paste(samples, "test of the median")
    } else if (measure == "bowley"){
      if(!is.null(p)){
        if(p > 0 & p < 0.5) u <- c(p, 0.5, 1 - p)
        else if(p > 0.5 & p < 1) u <- c(1 - p, 0.5, p)
        else stop("Argument p must be a numeric value in (0, 1) except 0.5.\n")
      } else u <- c(0.25, 0.5, 0.75)
      coef <- matrix(c(1, -2, 1, -1, 0, 1), nrow = 2, ncol = 3, byrow = TRUE)
      measure.name <- "Bowley's skew"
      method <- paste(samples, "test of Bowley's quantile skew measure")
    } else if (measure == "kelly"){
      u <- c(0.1, 0.5, 0.9)
      coef <- matrix(c(1, -2, 1, -1, 0, 1), nrow = 2, ncol = 3, byrow = TRUE)
      measure.name <- "Kelly's skew"
      method <- paste(samples, "test of Kelly's quantile skew measure")
    } else if (substr(measure, start = 1, stop = 5) == "groen"){
      if ((g.skew <- substr(measure, start = 6, stop = 6)) %in% c("R", "L")){
      } else stop("Error with Groeneveld's skew.  Use one of 'groenR' or 'groenL'.\n")
      if(!is.null(p)){
        if(p > 0 & p < 0.5) u <- c(p, 0.5, 1 - p)
        else if(p > 0.5 & p < 1) u <- c(1 - p, 0.5, p)
        else stop("Argument p must be a numeric value in (0, 1) except 0.5.\n")
      } else u <- c(0.25, 0.5, 0.75)
      if(g.skew == "R") denom <- c(-1, 1, 0)
      else denom <- c(0, 1, -1)
      coef <- matrix(c(1, -2, 1, denom), nrow = 2, ncol = 3, byrow = TRUE)
      measure.name <- "Groeneveld and Meeden's skew"
      method <- paste(samples, "test of Groeneveld and Meeden's quantile skew measure")
    } else if (measure == "moors"){
      u <- c(1/8, 2/8, 3/8, 5/8, 6/8, 7/8)
      coef <- matrix(c(-1, 0, 1, -1, 0, 1, 0, -1, 0, 0, 1, 0), nrow = 2, ncol = 6, byrow = TRUE)
      measure.name <- "Moors kurtosis"
      method <- paste(samples, "test of Moors kurtosis")
    } else if (measure == "lqw"){
      if(!is.null(p)){
        if(p > 0 & p < 0.5) u <- c(p/2, 0.25, (1 - p)/2)
        else stop("Argument p must be a numeric value in (0, 1/2).\n")
      } else u <- c(0.25/2, 0.25, 0.75/2)
      coef <- matrix(c(1, -2, 1, -1, 0, 1), nrow = 2, ncol = 3, byrow = TRUE)
      measure.name <- "Left quantile tail weight"
      method <- paste(samples, "test of left quantile tail weight")
    } else if (measure == "rqw"){
      if(!is.null(p)){
        if(p > 0.5 & p < 1) u <- c(1 - p/2, 0.75, (1 + p)/2)
        else stop("Argument p must be a numeric value in (1/2, 1).\n")
      } else u <- c(1 - 0.75/2, 0.75, (1 + 0.75)/2)
      coef <- matrix(c(1, -2, 1, -1, 0, 1), nrow = 2, ncol = 3, byrow = TRUE)
      measure.name <- "Right quantile tail weight"
      method <- paste(samples, "test of right quantile tail weight")
    } else if (substr(measure, start = 1, stop = 2) == "qr"){
      num <- substring(measure, 3)
      ratio.error <- "For quantile ratios, measure must be in format qrxxyy where xx and yy are integer numbers.\n"
      if(nchar(num) != 4) stop(ratio.error)
      coef <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = TRUE)
      u <- as.numeric(substring(measure, c(3, 5), c(4, 6)))/100
      if(any(is.na(u))) stop(ratio.error)
      measure.name <- paste0("ratio of ", u[1], " and ", u[2], " quantiles")
      method <- paste(samples, "test of ratio of quantiles")
    } else stop("Unknown choice for measure.\n")
  }
  else if (is.vector(coef)){
    if(!is.numeric(coef))
      stop("If coef is a vector, then it needs to be a numeric vector.\n")
    if(is.null(u))
      stop("Argument u required if coef is a numeric vector.\n")
    if (length(u) != length(coef))
      stop("Length of u needs to be equal to the length of coef.\n")
    if(length(coef) > 1) {
      measure.name <- "LCQ"
      method <- paste(samples, "test of a linear combination of quantiles (LCQ)")
    } else{
      measure.name <- paste0("quantile(u=", u,")")
      method <- paste(samples, "test of a single", measure.name)
    }
  }
  else if (is.matrix(coef)){
    if(!is.numeric(coef))
      stop("If coef is a matrix, then it needs to be a numeric matrix.\n")
    if(is.null(u))
      stop("Argument u required if coef is a numeric matrix.\n")
    if (length(u) != ncol(coef) | nrow(coef) != 2)
      stop("Matrix coef needs to have dimensions: ncol(coef)=length(u) and nrow(coef)=2.\n")
    measure.name <- "Ratio of LCQs"
    method <- paste(samples, "test of a ratio of two linear combinations of quantiles (LCQs)")
  }
  else stop("Argument coef must be either a character string, or a numeric vector or matrix.\n")

  if(!log.transf & is.matrix(coef)){
    if(nrow(coef) == 2) warning("You may wish to consider using a log transformation for ratios
                                (e.g. log.transf = TRUE).  If you choose to use
                                a log transformation you can also back-transform to the ratio
                                scale using (back.transf = TRUE) if you wish.\n")
  }

  alpha <- 1 - conf.level
  crit <- qnorm(1 - alpha/2)
  coef.string <- is.character(coef)

  if(! var.method %in% c("qor", "density"))
    stop("Argument method must be either 'qor' or 'density'.\n")

  qestx <- quantile(x, u, type = quantile.type)
  covQx <- qcov(x, u, method = var.method, quantile.type = quantile.type)$cov

  if (is.vector(coef)){
    estx <- sum(coef*qestx)
    sterrx <- sqrt(t(coef)%*%covQx%*%coef)[1, 1]
  } else{
    u1qx <- sum(coef[1, ]*qestx)
    u2qx <- sum(coef[2, ]*qestx)
    s1x <- (t(coef[1, ])%*% covQx %*%coef[1, ])[1, 1]
    s2x <- (t(coef[2, ])%*% covQx %*%coef[2, ])[1, 1]
    s12x <- (t(coef[1, ])%*% covQx %*%coef[2, ])[1, 1]

    estx <- u1qx/u2qx
    sterrx <- sqrt(s1x/u2qx^2 + u1qx^2/u2qx^4*s2x - 2*u1qx*s12x/u2qx^3)
  }
  if (log.transf){
    sterrx <- sterrx/estx
    estx <- log(estx)
    if (back.transf) transf.text <- NULL
    else  transf.text <- "(log transformed)"
  } else transf.text <- NULL

  if(!is.null(y)){
    qesty <- quantile(y, u, type = quantile.type)
    covQy <- qcov(y, u, method = var.method, quantile.type = quantile.type)$cov

    if (is.vector(coef)){
      esty <- sum(coef*qesty)
      sterry <- sqrt(t(coef)%*%covQy%*%coef)[1, 1]
    } else{
      u1qy <- sum(coef[1, ]*qesty)
      u2qy <- sum(coef[2, ]*qesty)
      s1y <- (t(coef[1, ])%*% covQy %*%coef[1, ])[1, 1]
      s2y <- (t(coef[2, ])%*% covQy %*%coef[2, ])[1, 1]
      s12y <- (t(coef[1, ])%*% covQy %*%coef[2, ])[1, 1]

      esty <- u1qy/u2qy
      sterry <- sqrt(s1y/u2qy^2 + u1qy^2/u2qy^4*s2y - 2*u1qy*s12y/u2qy^3)
    }


    if (log.transf){
      sterry <- sterry/esty
      esty <- log(esty)
      if (back.transf) transf.text <- NULL
      else  transf.text <- "(log transformed)"
    } else transf.text <- NULL

  }

  if (samples == "Two sample"){
    est <- estx - esty
    sterr <- sqrt(sterrx^2 + sterry^2)
  } else{
    est <- estx
    sterr <- sterrx
  }

  test.stat <- (est - true.q)/sterr
  names(test.stat) <- "Z"

  if(alternative == "less"){
    pval <- pnorm(test.stat)
    ci <- c(min.q, est + qnorm(conf.level)*sterr)
  } else if (alternative == "greater"){
    pval <- pnorm(test.stat, lower.tail = FALSE)
    ci <- c(est - qnorm(conf.level)*sterr, Inf)
  } else{
    pval <- 2*(1 - pnorm(abs(test.stat)))
    ci <- est + c(-1, 1)*crit*sterr
  }

  if(log.transf & back.transf){
    est <- exp(est)
    ci <- exp(ci)
    if(true.q == 0) true.q <- 1
  }

  attr(ci, 'conf.level') <- conf.level


  names(est) <- measure.name
  if(samples == "Two sample"){
    if (names(est) != "Ratio of LCQs") names(est) <- paste0(names(est), "s")
    if(log.transf & back.transf) names(est) <- paste("ratio of", names(est))
    else names(est) <- paste("difference in", names(est))
  }

  names(est) <- paste(names(est), transf.text)
  names(true.q) <- names(est)

  qres <- list(method = method,
               data.name = dname,
               statistic = test.stat,
               parameter = NULL,
               p.value = pval,
               alternative = alternative,
               estimate = est,
               null.value = true.q,
               conf.int = ci)

  class(qres) <- "htest"
  return(qres)
}
