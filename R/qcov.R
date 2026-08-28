#' Approximate Covariance Matrix Estimation for Vectors of Quantile Estimators
#' @description
#' compute a covariance matrix consisting of variances (on the diagonal) for quantile estimates and covariances (off-diagonal) between different quantile estimates
#' @details
#'
#' This function computes an approximate covariance matrix for a vector of sample quantile estimators. The covariance structure is based on the asymptotic relationship, for \eqn{u_i\leq u_j} and sample size \eqn{n},
#' \deqn{\mathrm{Cov}(\widehat Q(u_i), \widehat Q(u_j)) \approx
#' \displaystyle\frac{u_i(1 - u_j)}{n}\, q(u_i)\, q(u_j),}
#' where \eqn{q(u)=Q'(u)} is the quantile density.
#'
#' The quantile density is estimated by \code{qden}. By default, \code{method = "qor"} is used, which estimates the quantile density using the quantile optimality ratio (QOR) approach of Prendergast and Staudte (2016). Alternatively, \code{method = "density"} estimates the quantile density by fitting a density function and taking its reciprocal at the estimated quantiles.
#'
#' The argument \code{dist} specifies the working distribution used by \code{qden}. It may be one of the built-in distributions supported by \code{qden}, or a quantile function name when a user-supplied distribution is to be used. If \code{params} is \code{NULL}, the required distributional parameters are estimated from the data where appropriate.
#'
#' If \code{method = "qor"}, the bandwidth is selected using the QOR approach. The argument \code{bw.correct} controls the boundary correction used in that bandwidth calculation.
#'
#' If \code{method = "density"}, additional arguments in \code{...} are passed to \code{density} when density estimation is used directly.
#' @param x a numeric vector of data values.
#' @param u a numeric vector of probability values in the interval \eqn{(0,1)} specifying the quantiles to be estimated. Missing values are not allowed.
#' @param method the approach used to estimate the quantile density function. Either \code{"qor"} or \code{"density"}.
#' @param dist a character string naming a supported distribution or a quantile function name used by \code{qden}.
#' @param quantile.type argument for the quantile function. Default is set to \code{8} so that output is consistent with default quantile function use and other functions such as \code{IQR} (see \code{?quantile} for details).
#' @param bw.correct logical; if \code{TRUE}, the bandwidth is corrected near the boundary.
#' @param params a list of parameter values for the chosen distribution. If \code{NULL}, parameters are estimated where needed.
#' @param ... additional arguments to be passed to \code{qden} when \code{method = "qor"} and to \code{density} when \code{method = "density"}.
#' @return a covariance matrix consisting of variances (on the diagonal) for quantile estimates and covariances (off-diagonal) between different quantile estimates.
#' @references
#' Prendergast, L. A., & Staudte, R. G. (2016). Exploiting the quantile optimality ratio in finding confidence intervals for quantiles. \emph{Stat}, \strong{5}(1), 70--81.
#'
#' Prendergast, L. A., Dedduwakumara, D. S., & Staudte, R. G. (2024). \emph{rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles}. Preprint, pages 1--13.
#' @seealso
#' \code{\link{qden}} for estimating the quantile density function,
#' \code{\link{qor}} for quantile optimality ratio values, and
#' \code{\link{qrcov}} for covariance matrices of ratios of linear combinations of quantiles.
#' @export
#' @examples
#' set.seed(1234)
#' x <- rnorm(100)
#'
#' # Covariance matrix for sample quartiles
#' qcov(x, c(0.25, 0.5, 0.75))
#'
#' # Density-based version
#' qcov(x, c(0.25, 0.5, 0.75), method = "density")

qcov <- function(x, u,
                 method = "qor",
                 dist = "gl",
                 quantile.type = 8,
                 bw.correct = TRUE,
                 params = NULL,
                 ...)
{
  if (!is.numeric(x))
    stop("Argument 'x' must be numeric.")

  if(any(u <= 0 | u >=1) | anyNA(u)){
    stop("Argument u must be a numeric vector of probability values between, but not including, 0 and 1.")
  }

  n <- length(x)

  u1u <- u %*% t(1 - u)
  u1u <- pmin(u1u, t(u1u))

  if (method == "qor") {

    qden.hat <- qden(x, u,dist = dist,method = "qor",bw.correct = bw.correct,params = params,...)
    covQ <- u1u * tcrossprod(qden.hat) / n

  } else if (method == "density") {

    qden.hat <- qden(x, u,dist = dist,method = "density",quantile.type = quantile.type,params = params,...)
    covQ <- u1u * tcrossprod(qden.hat) / n

  } else {
    stop("'method' must be either 'qor' or 'density'.\n")
  }

  rownames(covQ) <- u
  colnames(covQ) <- u

  return(covQ)
}
