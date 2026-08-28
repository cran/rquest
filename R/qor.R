#' Quantile Optimality Ratio
#' @description
#' compute the quantile optimality ratio (QOR) for a vector of probabilities from either a named distribution or a user-supplied quantile function
#' @details
#' This function evaluates the quantile optimality ratio (QOR, Prendergast & Staudte, 2016) for probabilities in \code{u}.
#' For a quantile function \eqn{Q}, let \eqn{q(u)=Q'(u)} denote the quantile density. The QOR is
#' \deqn{\mathrm{QOR}(u)=\frac{q(u)}{q''(u)}=\frac{Q'(u)}{Q'''(u)}.}
#' If \code{dist} is a character string, the QOR is computed using a closed-form expression whenever one is available.
#' For \code{"norm"} and \code{"exp"}, the QOR does not depend on the distribution parameters and therefore \code{x}
#' is not required. For \code{"lnorm"}, the required parameter is estimated from \code{x} when \code{params} is \code{NULL}.
#' For \code{"gl"}, the FKML generalized lambda distribution is fitted using \code{fit.fkml}
#' and the options supplied in \code{gl.control}, unless \code{params}
#' is supplied directly.
#'
#' If \code{dist} is a function, it is treated as a quantile function. In that case, \code{params} supplies the arguments
#' passed to the function. If \code{params} is \code{NULL}, the default parameter values of the supplied quantile function
#' are used and a warning is issued. For user-supplied quantile functions, the QOR is evaluated numerically using finite
#' differences.
#' @param u a numeric vector of probability values in the interval \eqn{(0,1)}. Missing values are not allowed.
#' @param dist a character string naming a supported distribution, or a function giving a quantile function.
#' Supported character values currently include \code{"gl"}, \code{"lnorm"}, \code{"norm"}, and \code{"exp"}.
#' @param x a numeric vector of data values. This is used to estimate distribution parameters when \code{params} is \code{NULL}
#' for supported built-in distributions.
#' @param params a list of parameter values to use. If \code{NULL} and \code{dist} is a supported built-in distribution,
#' the parameters are estimated from \code{x} where needed. If \code{NULL} and \code{dist} is a function, the default
#' arguments of that quantile function are used with a warning.
#' @param gl.control a list of control arguments passed to \code{fit.fkml} when \code{dist = "gl"}.
#' @return a list with components \code{qor} (the QOR values at \code{u}) and \code{params} (the parameter values used).
#' @seealso
#' \code{\link{qcov}} for estimating the covariance matrix of sample quantiles,
#' \code{\link{qrcov}} for estimating covariance matrices of ratios of linear
#' combinations of quantiles, and \code{\link[gld]{fit.fkml}} from package
#' \pkg{gld} for fitting the FKML generalized lambda distribution.
#'
#' @references
#' Prendergast, L. A., & Staudte, R. G. (2016). Exploiting the quantile optimality ratio in finding confidence intervals for quantiles. Stat, 5(1), 70-81
#'
#' @export
#' @examples
#' set.seed(1234)
#' x <- rlnorm(100)
#'
#' # QOR using the flexible FKML GLD distribution (default)
#' qor(seq(0.1, 0.9, by = 0.2), x = x)
#'
#' # QOR for the lognormal distribution
#' qor(seq(0.1, 0.9, by = 0.2), dist = "lnorm", x = x)
#'
#' # QOR for the exponential distribution
#' qor(seq(0.1, 0.9, by = 0.2), dist = "exp")
#'
#' # QOR for a user-supplied quantile function
#' qor(c(0.2, 0.5, 0.8),
#'     dist = qbeta,
#'     params = list(shape1 = 2, shape2 = 5))
qor <- function(u, dist = "gl", x = NULL, params = NULL, gl.control = list(method = "Lmom")) {
  u <- as.numeric(u)
  if (anyNA(u) || any(u <= 0 | u >= 1)) {
    stop("Argument 'u' must contain probability values between, but not including, 0 and 1.")
  }
  out <- NULL

  if (is.character(dist)) {
    if (dist %in% c("gl", "lnorm") && is.null(x) && is.null(params))
      stop("Must provide 'x' if 'params' is NULL.\n")

    if(dist == "gl"){
      if (is.null(params)) {
        res <- do.call(gld::fit.fkml, c(list(x = x), gl.control))
        l <- res$lambda
        params <- as.list(l)
        names(params) <- c("lambda1", "lambda2", "lambda3", "lambda4")
      } else {
        l <- unlist(params)
      }
      num <- u^(l[3] - 1) + (1 - u)^(l[4] - 1)
      den <- (l[3] - 1)*(l[3] - 2)*u^(l[3] - 3) + (l[4] - 2)*(l[4] - 1)*(1 - u)^(l[4] - 3)
      out <- num/den
    } else if (dist == "lnorm") {
      if (!is.null(x) && any(x <= 0)) {stop("Values in 'x' must be positive if using 'lnorm'.\n")}
      if (is.null(params)) {
        sdlog <- sd(log(x))
        params <- list(sdlog = sdlog)
      } else {
        sdlog <- params$sdlog
      }

      zu <- qnorm(u)
      out <- dnorm(zu)^2 /
        (1 + sdlog^2 + 3 * sdlog * zu + 2 * zu^2)
    } else if (dist == "norm"){
      zu <- qnorm(u)
      out <- dnorm(zu)^2 / (1 + 2 * zu^2)
    } else if (dist == "exp"){
      out <- (1 - u)^2/2
    } else {
        stop("Unknown argument for 'dist'.\n")
    }
  } else if (is.function(dist)) {

    qf <- match.fun(dist)
    Q <- function(u) do.call(qf, c(list(p = u), params))

    if (is.null(params)) {
      warning("Using the default parameter values of the supplied quantile function.\n")
      params <- list()
    }

    h <- pmin(0.0001, u/4, (1 - u)/4)
    d1 <- (Q(u + h) - Q(u - h))/(2 * h)
    d3 <- (Q(u + 2*h) - 2*Q(u + h) + 2*Q(u - h) - Q(u - 2*h))/(2*h^3)
    out <- d1/d3
  } else {
    stop("'dist' must be a character string or a function.\n")
  }
  list(qor = out, params = params)
}
