#' qcov
#' @description
#' compute a covariance matrix consisting of variances (on the diagonal) for quantile estimates and covariances (off-diagonal) between different quantile estimates
#' @details
#' This function computes a covariance matrix for a vector of quantile estimators.
#' This is done via estimating the inverted density function evaluated at the respective quantiles.
#' The default for this is to use the quantile optimality ratio (QOR) approach (Prendergast & Staudte, 2016) which computes an optimal bandwidth.
#' Alternatively, using `method = "density"` will use the generic density function to estimate the density.
#' @param x a numeric vector of data values.
#' @param u a numeric vector of probability values in \[0, 1\]) specifying the quantiles to be estimate.
#' @param method approach use to estimate the quantile density function. Either "qor" or "density".
#' @param FUN QOR function for the log-normal
#' @param quantile.type argument for the quantile function.  Default is set to 8 so that output is consistent with default quantile function use and other functions such as IQR (see help file for `quantile()`
#' for more details)
#' @param bw.correct replace bw by the values of v when v<=bw (see Prendergast & Staudte (2016b) for more details)
#' @return a covariance matrix consisting of variances (on the diagonal) for quantile estimates and covariances (off-diagonal) between different quantile estimates
#' @references
#' Prendergast, L. A., & Staudte, R. G. (2016). Exploiting the quantile optimality ratio in finding confidence intervals for quantiles. Stat, 5(1), 70-81
#'
#' Prendergast, L. A., Dedduwakumara, D.S. & Staudte, R.G. (2024) rquest: An R package for hypothesis tests and confidence intervals
#' for quantiles and summary measures based on quantiles, preprint, pages 1-13
#'
#' @export
#'
#' @examples
#' # Create some data
#' set.seed(1234)
#' x <- rnorm(100)
#'
#' # Compute the variance-covariance matrix for sample quartiles.
#' qcov(x, c(0.25, 0.5, 0.75))

qcov <- function(x, u,
                 method = "qor",
                 FUN = qor.ln,
                 quantile.type = 8,
                 bw.correct = TRUE){

  n <- length(x)
  qest <- quantile(x, u, type = quantile.type)

  u1u <- u %*% t(1 - u)
  u1u <- pmin(u1u, t(u1u))

  if (method == "qor"){
    qor <- FUN(u)
    bw <- 15^(1/5)*abs(qor)^(2/5)/n^(1/5)
    if (bw.correct) bw[u <= bw] <- u[u <= bw]#replace bw by the values of v when v<=bw

    kernepach <- function(u) 3/4*(1 - u^2)*(abs(u) <= 1)
    J <- length(u)
    m1 <- matrix(u, nrow = J, ncol = n, byrow = FALSE)
    m2 <- matrix(1:n, nrow = J, ncol = n, byrow = TRUE)

    consts <- kernepach((m1 - (m2 - 1)/n)*(1/bw))*(1/bw) -
      kernepach((m1 - m2/n)*(1/bw))*(1/bw)

    x.sorted <- sort(x)
    q.hat <- c(consts%*%x.sorted)

    covQ <- u1u*(q.hat %*% t(q.hat))/n

  } else if (method == "density"){
    dest <- density(x)
    df <- approxfun(dest)
    covQ <- u1u*((1/df(qest)) %*% t(1/df(qest)))/n
  }

  rownames(covQ) <- u
  colnames(covQ) <- u

  list(cov = covQ)
}
