#' Hypothesis Tests and Confidence Intervals for Quantile-based Inequality Measures
#' @description
#' carry out hypothesis tests and obtain associated confidence intervals for quantile based inequality measures
#' @details
#'
#' This function performs hypothesis testing and calculates the corresponding confidence intervals for inequality measures based on quantiles.
#' Let \eqn{l_1(p)} and \eqn{l_2(p)} be functions of \eqn{p \in [0,1]} that defined linear combinations of quantiles.  For example, if
#' \eqn{x_{0.4}} and \eqn{x_{0.6}} are the 0.4 and 0.6 quantiles respectively, and the ratio \eqn{x_{0.4}/(x_{0.4} + x_{0.6})} is of
#' interest, then possible choices are \eqn{l_1(p)=x_p} and \eqn{l_2(p)=x_p + x_{1-p}} for \eqn{p=0.4}.  The inequality measures supported
#' are of the form
#' \deqn{I = 1 - \displaystyle\int^1_0\frac{l_1(p)}{l_2(p)}dp}
#' and where \eqn{I\in [0, 1]}.
#'
#' Estimation of \eqn{I} is done numerically over a grid of points determined by integer \eqn{J} (default 100), where \eqn{p_j = (j-1/2)/J},
#' and the estimate to \eqn{I} is
#' \deqn{\widehat{I} = 1 - \frac{1}{J}\displaystyle\sum_{j=1}^J\frac{\hat{l}_1(p_j)}{\hat{l}_2(p_j)}}
#' where \eqn{\hat{l}_1(p_j)} and \eqn{\hat{l}_2(p_j)} are the estimated linear combinations of quantiles at \eqn{p_j}.
#'
#' The available options for quantile based measures in argument `measure` are shown below.
#'
#' * `"QRI"`                : Quantile Ratio Index from Prendergast & Staudte (2018). This is the default choice.
#' * `"G1"`, `"G2"`, `"G3`" :  Quantile variants of the Gini index from Prendergast & Staudte (2016a).
#' * `"S1"`, `"S2"`         :  'Strategy 1' and 'Strategy 2' income inequalities from Brazauskas et al. (2024).
#'
#' The default `var.method="qor"` is to estimate the probability density function directly using the lognormal Quantile Optimality Ratio (QOR)
#' for choosing a suitable bandwidth (Prendergast & Staudte,2016b). Alternatively, the variances can be
#' estimated by inverting a density estimator evaluated at the quantiles and this can be done using `var.method = "density"`. If `var.method = "density"`,
#' then the function density is used to estimate the probability density function which is needed for the calculation of the covariance matrix using function qcov.
#' If needed, additional arguments can be passed to density (see ?density for details on possible additional arguments).
#'
#' It is also possible for users to define their own inequality measures either using a formula or by creating a list.  The formula approach is
#' the most straightforward of the two and can be defined in terms `Q` to specify a quantile function, `p` and then specifying the
#' formula which must be a ratio of linear combinations of quantiles, and where the coefficients may also contain `p`.  As an example,
#' the QRI measure would be `measure = ~ Q(p/2)/Q(1 - p/2)` while G2 would require `measure = ~ 2*Q(p/2)/Q(1 - p/2)`.  More details
#' and examples can be found in Prendergast et al. (2024).
#'
#' If `measure` is a list, it must be of the form `list(u=...,coef1=...,coef2)` where `u` is a numeric vector indicating the
#' probabilities for every quantile estimate needed for \eqn{\hat{I}}, `coef1` and `coef2` are numeric matrices with `J` rows
#' and number of columns equal to the length of `u`.  Each row of `coef1` stores the coefficients for the numerator
#' linear combination, and ditto `coef2` for the denominator.
#'
#' For more information and further examples, including user defined measures, see Prendergast, Dedduwakumara & Staudte (2024).
#'
#' @param x a numeric vector of data values.
#' @param y an optional second vector of data values for two-sample testing.
#' @param J number of grid points
#' @param measure character string of quantile based inequality measure to be estimated or a list of parameter choices or formula for user defined measures (see details for available measures and user defined measures).
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param quantile.type argument for the quantile function.  Default is set to 8 so that output is consistent with default quantile function use and other functions such as IQR (see help file for `quantile()`
#' for more details)
#' @param var.method approach use to estimate the quantile density function.  Either "qor"(default) or "density".(See details).
#' @param conf.level coverage for the estimated confidence interval.
#' @param true.ineq the specified hypothesized value of the inequality measure or the difference of the inequality depending on whether it was a one-sample test or a two-sample test.
#' @param ... additional arguments to be passed to function qcov when var.method = “density” is used.
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
#' Brazauskas, V., Greselin, F., & Zitikis, R. (2024). Measuring income inequality via percentile relativities. Quality & Quantity, 58(5), 4859-4896.
#'
#' @export
#'
#' @examples
#' # Create some data
#' x <- c(8.43,7.08,8.79,8.88,7.87,5.94,8.79,5.46,8.11,7.08)
#' y <- c(13.44,13.65,14.77,9.51,14.07,10.92,11.59,13.42,8.93,10.88)
#'
#' # Test of equality of QRI measure between two groups
#' qineq(x, y, measure = "QRI")
#'
#' # Another way by defining the formula for the ratio in QRI
#' qineq(x, y, measure = ~ Q(p/2)/Q(1 - p/2))
#'
#' # Another way by defining the quantiles and coefficients directly in a list
#' J <- 100
#' p <- (1:J - 0.5)/J
#' u <- sort(c(p/2, 1 - p/2))
#' num <- cbind(diag(rep(1, J)), matrix(0, J, J))
#' den <- cbind(matrix(0, J, J), diag(J)[, J:1])
#' qineq(x, y, measure = list(u = u, coef1 = num, coef2 = den))
#'

qineq <- function (x, y = NULL, J = 100, measure = "QRI", alternative = c("two.sided",
                                                                          "less", "greater"), quantile.type = 8, var.method = "qor",
                   conf.level = 0.95, true.ineq = 0.5, ...)
{
  if (!is.numeric(x))
    stop("Argument 'x' must be numeric.")

  alternative <- match.arg(alternative)
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

  measures <- c("QRI", "G1", "G2", "G3", "S1", "S2")
  p <- (1:J - 0.5)/J
  if(is.character(measure)){
    if (!is.na(match(measure, measures))) {
      params <- switch(measure,
                       "QRI" = list(u      = sort(c(p/2, 1 - p/2)),
                                    coef1  = cbind(diag(rep(1, J)), matrix(0, J, J)),
                                    coef2  = cbind(matrix(0, J, J), diag(J)[, J:1])),
                       "G1"  = list(u      = c(p/2, 0.5),
                                    coef1  = cbind(2*diag(p), rep(0, J)),
                                    coef2  = cbind(matrix(0, J, J), rep(1, J))),
                       "G2"  = list(u      = sort(c(p/2, 1 - p/2)),
                                    coef1  = cbind(2*diag(p), matrix(0, J, J)),
                                    coef2  = cbind(matrix(0, J, J), diag(J)[, J:1])),
                       "G3"  = list(u      = sort(c(p/2, 1 - p/2)),
                                    coef1  = cbind(4*diag(p), matrix(0, J, J)),
                                    coef2  = cbind(diag(J), diag(J)[, J:1])),
                       "S1"  = list(u      = c(p/2, 0.5),
                                    coef1  = cbind(diag(J), rep(0, J)),
                                    coef2  = cbind(matrix(0, J, J), rep(1, J))),
                       "S2"  = list(u      = c(p/2, 0.5 + p/2),
                                    coef1  = cbind(diag(J), matrix(0, J, J)),
                                    coef2  = cbind(matrix(0, J, J), diag(J))))
    } else stop(
      "Unknown inequality measure specified.  Please ensure 'measure' is one of: ",
      paste(sprintf("'%s'", measures), collapse = ", ")
    )
  } else if (rlang::is_formula(measure)){
    params <- get.coefs(measure, p)
    "measure" <- "user defined"
  } else if (is.list(measure)) {

    if (!all(c("u", "coef1", "coef2") %in% names(measure)))
      stop("If you are defining your own inequality measure, then 'measure' must be a list containing elements 'u', 'coef1' and 'coef2'.")

    u <- measure$u
    if (!is.numeric(u))
      stop("'u' in list 'measure' must be a numeric.")
    if(any(u <= 0 | u >=1) | anyNA(u)){
      stop("'u' in list 'measure' must be a numeric vector of probability values between, but not including, 0 and 1.")
    }
    coef1 <- measure$coef1
    coef2 <- measure$coef2
    if (!is.matrix(coef1) | !is.matrix(coef2))
      stop("Both 'coef1' and 'coef2' in list 'measure' must be a matrix.")

    len.u <- length(u)
    if(!all(dim(coef1) == dim(coef2)) | ncol(coef1) != len.u | ncol(coef2) != len.u){
      stop("Both 'coef1' and 'coef2' in list 'measure' must be numeric matrices with the same dimension, and the number of columns equal to the length of 'u'.")
    }
    params <- measure
    "measure" <- "user defined"
  } else{
    stop("Argument 'measure' must be either character string or a list.")
  }

  measure <- paste0(measure, " statistic")
  alpha <- 1 - conf.level
  crit <- qnorm(1 - alpha/2)
  u <- params$u
  coef1 <- unname(params$coef1)
  coef2 <- unname(params$coef2)

  qrcov.resultx <- qrcov(x, u = u, coef1 = coef1, coef2 = coef2, quantile.type = quantile.type, method = var.method, ...)
  J <- length(qrcov.resultx$ratios)
  estx <- 1 - sum(qrcov.resultx$ratios)/J
  sterrx <- sqrt(sum(qrcov.resultx$cov)/J^2)
  names(estx) <- measure
  method <- paste(samples, "test of the", measure)

  if (!is.null(y)) {
    qrcov.resulty <- qrcov(y, u = u, coef1 = coef1, coef2 = coef2, quantile.type = quantile.type,method = var.method, ...)
    esty <- 1 - sum(qrcov.resulty$ratios)/J
    sterry <- sqrt(sum(qrcov.resulty$cov)/J^2)
    names(esty) <- measure
  }


  if (samples == "Two sample") {
    est <- estx - esty
    sterr <- sqrt(sterrx^2 + sterry^2)
    names(est) <- paste0("difference in ", measure)
    if (true.ineq == 0.5)
      true.ineq <- 0
  }
  else {
    est <- estx
    sterr <- sterrx
  }
  test.stat <- (est - true.ineq)/sterr
  names(test.stat) <- "Z"
  if (alternative == "less") {
    pval <- pnorm(test.stat)
    if (samples == "Two sample") {
      ci <- c(-1, est + qnorm(conf.level) * sterr)
    } else {
      ci <- c(0, est + qnorm(conf.level) * sterr)
    }
  }
  else if (alternative == "greater") {
    pval <- pnorm(test.stat, lower.tail = FALSE)
    ci <- c(est - qnorm(conf.level) * sterr, 1)
  }
  else {
    pval <- 2 * (1 - pnorm(abs(test.stat)))
    ci <- est + c(-1, 1) * crit * sterr
  }
  attr(ci, "conf.level") <- conf.level
  names(true.ineq) <- names(est)
  ineqres <- list(method = method, data.name = dname, statistic = test.stat,
                  parameter = NULL, p.value = pval, alternative = alternative,
                  estimate = est, null.value = true.ineq, conf.int = ci)
  class(ineqres) <- "htest"
  return(ineqres)
}
