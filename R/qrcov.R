#' Approximate Covariance Matrix Estimation for Ratios of Linear Combinations of Quantile Estimators
#' @description
#' compute a covariance matrix consisting of variances (on the diagonal) for ratios of linear combinations of quantile estimates
#' and covariances (off-diagonal) between different ratios
#' @details
#' This function computes a sample covariance matrix for ratios of linear combinations of quantile estimates from a single sample.
#' This is done via a sample covariance matrix for the quantiles used in the linear combinations which is obtained using function `qcov`.
#' Let \eqn{\mathbf{q}} denote a vector of quantile estimators, \eqn{\mathbf{S}} denote the covariance matrix for \eqn{\mathbf{q}},
#' and \eqn{\mathbf{c}_{1i}} and \eqn{\mathbf{c}_{2i}} \eqn{i=1,\ldots,r} be the numerator and denominator vectors of coefficients
#' that define the linear combinations used in each of the \eqn{r} ratios (e.g. the \eqn{i}th ratio is
#' \eqn{R_i = (\mathbf{c}_{1i}^\top \mathbf{q})/(\mathbf{c}_{2i}^\top \mathbf{q})}).  Then, using the Delta method, an approximate
#' covariance matrix for \eqn{[R_1,\ldots,R_r]^\top} is
#' \deqn{\mathbf{J}\mathbf{S}\mathbf{J}^\top} where \eqn{\mathbf{J}} is the Jacobian matrix whose \eqn{i}th row is the gradient
#' vector of \eqn{R_i}.

#' @param x a numeric vector of data values.
#' @param u a numeric vector of probability values in the interval (0,1) specifying the quantiles to be estimated. Note that u must include numeric values between, and not including, 0 and 1 and missing values are not allowed.
#' @param coef1 a numeric matrix whose \eqn{i}th row contains the coefficients for the numerator linear combination for the \eqn{i}th ratio.
#' @param coef2 a numeric matrix whose \eqn{i}th row contains the coefficients for the denominator linear combination for the \eqn{i}th ratio.
#' @param quantile.type argument for the quantile function.  Default is set to 8 so that output is consistent with default quantile function use and other functions such as IQR (see help file for `quantile()`
#' for more details)
#' @param ... additional arguments to be passed to function density when method = “density” is used.
#' @return a list consisting of `ratios` (the ratio estimates) and `cov` (the estimated covariance matrix).
#' @references
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
#' # Compute the variance-covariance matrix for the ratio of third
#' # quartile and the median and the ratio of first quartile and the median.
#' coef1 <- matrix(c(0, 0, 1, 1, 0, 0), nrow = 2, byrow = TRUE)
#' coef2 <- matrix(c(0, 1, 0, 0, 1, 0), nrow = 2, byrow = TRUE)
#' qrcov(x, c(0.25, 0.5, 0.75), coef1 = coef1, coef2 = coef2)


qrcov <- function(x, u, coef1, coef2, quantile.type = 8, ...) {
  if (!is.numeric(x))
    stop("Argument 'x' must be numeric.")
  if (any(u <= 0 | u >= 1) | anyNA(u)) {
    stop("Argument 'u' must be a numeric vector of probability values between, but not including, 0 and 1.")
  }
  if(!is.numeric(coef1) | !is.numeric(coef2) | !is.matrix(coef1) | !is.matrix(coef2)){
    stop("Arguments 'coef1' and 'coef2' both must be a numeric matrix.")
  }
  len.u <- length(u)
  if(!all(dim(coef1) == dim(coef2)) | ncol(coef1) != len.u | ncol(coef2) != len.u){
    stop("Arguments 'coef1' and 'coef2' both must be numeric matrices with the same dimension, and the number of columns equal to the length of argument 'u'.")
  }
  qest <- quantile(x, u, type = quantile.type)
  covQ <- qcov(x, u, quantile.type = quantile.type, ...)
  p <- length(qest)
  m <- nrow(coef1)
  Jac <- matrix(0, nrow = m, ncol = p)
  ratios <- numeric(m)

  for (r in seq_len(m)) {
    c1r <- as.vector(coef1[r, ])
    c2r <- as.vector(coef2[r, ])
    num <- drop(crossprod(c1r, qest))
    den <- drop(crossprod(c2r, qest))
    ratios[r] <- num / den
    Jac[r, ] <- (c1r * den - c2r * num) / den^2
  }

  V <- Jac %*% covQ %*% t(Jac)
  dimnames(V) <- list(paste0("R", seq_len(m)), paste0("R", seq_len(m)))

  list(
    ratios = ratios,
    cov = V
  )
}
