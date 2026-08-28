#' Quantile Density Estimation
#' @description
#' estimate the quantile density function using either a direct estimation approach with a bandwidth selected by the quantile optimality ratio (QOR), or by the density method which uses the inverse of the estimated density at the nominated quantiles
#' @details
#' This function estimates the quantile density function at a vector of probability values \code{u}. Two estimation methods are available.
#'
#' If \code{method = "qor"}, a bandwidth is selected using the quantile optimality ratio (QOR) and a kernel-based estimator is used.
#' The bandwidth depends on the selected kernel through its variance and roughness constants. Supported kernels are \code{"epanechnikov"},
#' \code{"gaussian"}, \code{"rectangular"}, \code{"triangular"}, \code{"biweight"}, \code{"cosine"}, and \code{"optcosine"}.  By default, the QOR
#' method uses the flexible four parameter Generalized Lambda Distribution (GLD) to compute the QOR.  The method of L-moments is used for
#' GLD parameter estimation, although other estimators can be chosen by passing arguments through the \code{gl.control} argument to the \code{fit.fkml} function from the
#' \pkg{gld} package which is used for estimation.  Other included distributions for computing the QOR are the lognormal \code{dist = "lnorm"}, Gaussian \code{dist = "norm"},
#' and exponential \code{"exp"}.  If another character string is passed to \code{dist}, then the function will attempt to locate
#' the associated quantile function and so \code{dist} needs to follow R's naming convention of prefixing
#' the distribution name with q.  E.g., \code{dist = "beta"} will use the associated quantile function for the beta distribution, \code{qbeta}.  The parameter values to be used
#' for this distribution can be specified using the list argument \code{params}.  If \code{params} is \code{NULL}, the function
#' will attempt to use default choices for these parameters.
#'
#' If \code{method = "density"}, the quantile density is estimated as the reciprocal of a density estimate evaluated at the estimated quantiles.
#' When \code{params} is \code{NULL}, the density estimate is obtained directly from \code{density}. When \code{params} is supplied, the corresponding quantile and density functions from argument \code{dist} are used.
#' E.g., if \code{dist = "beta"}, then the quantile density will be computed using the \code{qbeta} and \code{dbeta} functions.
#'
#' @param x a numeric vector of data values.
#' @param u a numeric vector of probability values in the interval \eqn{(0,1)} specifying where the quantile density is to be estimated. Missing values are not allowed.
#' @param dist a character string naming a distribution or a quantile function name. Supported built-in distributions include \code{"gl"}, \code{"lnorm"}, \code{"norm"}, and \code{"exp"}.
#' @param method the estimation method to use. The options are \code{"qor"} and \code{"density"}.
#' @param kernel the kernel to use when \code{method = "qor"}.
#' @param params a list of parameter values. If \code{NULL}, parameters are estimated where needed.
#' @param gl.control a list of control arguments passed to \code{fit.fkml} when \code{dist = "gl"}.
#' @param quantile.type argument for the quantile function used when \code{method = "density"}. Default is \code{8} so that output is consistent with other functions such as \code{IQR}.
#' @param bw.correct logical; if \code{TRUE}, the bandwidth is corrected near the boundary.
#' @param large.n.switch logical; if \code{TRUE} and \code{dist = "gl"} with large \code{n}, then \code{"lnorm"} is used for bandwidth selection.
#' @param ... additional arguments to be passed to \code{density} when \code{method = "density"} is used, or additional arguments passed to \code{qor} when \code{method = "qor"} is used.
#' @return a vector of estimated quantile density values at \code{u}.
#' @seealso
#' \code{\link{qor}} for quantile optimality ratio values,
#' \code{\link{qcov}} for covariance estimation of sample quantiles, and \code{\link{qrcov}} for covariance estimation of ratios of linear combinations of quantiles.
#' @references
#' Prendergast, L. A., Dedduwakumara, D.S. & Staudte, R.G. (2024). \emph{rquest: An R package for hypothesis tests and confidence intervals for quantiles and summary measures based on quantiles}. Preprint, pages 1--13.
#'
#' @export
#' @examples
#' set.seed(1234)
#' x <- rnorm(100)
#'
#' # QOR-based quantile density estimation using the flexible GLD distribution (default method)
#' qden(x, c(0.25, 0.5, 0.75), method = "qor")
#'
#' # QOR-based quantile density estimation using the normal distribution
#' qden(x, c(0.25, 0.5, 0.75), dist = "norm", method = "qor")
#'
#' # Density-based quantile density estimation
#' qden(x, c(0.25, 0.5, 0.75), method = "density")
qden <- function(x, u, dist = "gl", method = "qor",
                       kernel = "epanechnikov",
                       params = NULL, gl.control = list(method = "Lmom"),
                       quantile.type = 8,
                       bw.correct = TRUE,
                       large.n.switch = TRUE,
                       ...){

  if (!is.numeric(x))
    stop("Argument 'x' must be numeric.")
  if(any(u <= 0 | u >=1) | anyNA(u)){
    stop("Argument u must be a numeric vector of probability values between, but not including, 0 and 1.")
  }

  n <- length(x)

  if(n > 1000 & dist == "gl" & large.n.switch){
    dist <- "lnorm"
    params <- NULL
    warning("For efficiency, the 'lnorm' distribution is being used for bandwidth selection.  To keep using 'gl' use 'large.n.switch = FALSE'.\n ")
  }
  if(method == "qor"){

    kernel <- match.arg(
      kernel,
      c("epanechnikov", "gaussian", "rectangular",
        "triangular", "biweight", "cosine", "optcosine")
    )

    Kinf <- switch(kernel,
                       "epanechnikov" = list(K = function(u) (3/4)*(1 - u^2)*(abs(u) <= 1),
                                             var = 1/5, rough = 3/5),
                       "gaussian"     = list(K = function(u) dnorm(u),
                                             var = 1, rough = 1/(2 * sqrt(pi))),
                       "rectangular"  = list(K = function(u) (1/2)*(abs(u) <= 1),
                                             var = 1/3, rough = 1/2),
                       "triangular"   = list(K = function(u) (1 - abs(u))*(abs(u) <= 1),
                                             var = 1/6, rough = 2/3),
                       "biweight"     = list(K = function(u) (15/16)*(1 - u^2)^2*(abs(u) <= 1),
                                             var = 1/7, rough = 5/7),
                       "cosine"       = list(K = function(u) (pi/4) * cos(pi*u/2) * (abs(u) <= 1),
                                             var = 1 - 8/pi^2, rough = pi^2/16),
                       "optcosine"    = list(K = function(u) (pi/4)*cos(pi*u/2)*(abs(u) <= 1),
                                             var = 1 - 8/pi^2, rough = pi^2/16)
    )
    if(! dist %in% c("gl", "lnorm", "norm", "exp")){
      dist <- get0(paste0("q", dist), mode = "function", inherits = TRUE)
    }
    qor.res <- qor(u, dist = dist, x = x, params = params, gl.control = gl.control, ...)$qor

    bw <- (Kinf$rough/Kinf$var^2)^(1/5) * abs(qor.res)^(2/5)/n^(1/5)
    if (bw.correct) bw[u <= bw] <- u[u <= bw]

    J <- length(u)
    m1 <- matrix(u, nrow = J, ncol = n, byrow = FALSE)
    m2 <- matrix(1:n, nrow = J, ncol = n, byrow = TRUE)

    consts <- Kinf$K((m1 - (m2 - 1)/n)*(1/bw))*(1/bw) -
      Kinf$K((m1 - m2/n)*(1/bw))*(1/bw)
    x.sorted <- sort(x)
    out <- c(consts %*% x.sorted)
  } else if (method == "density"){
    if(is.null(params)){
      dest <- density(x = x, ...)
      df <- approxfun(dest)
      qest <-  quantile(x, u, type = quantile.type)
      out <- 1/df(qest)
    } else{
      qfun <- get0(paste0("q", dist), mode = "function", inherits = TRUE)
      dfun <- get0(paste0("d", dist), mode = "function", inherits = TRUE)
      qpop <- do.call(qfun, c(list(p = u), params))
      out <- 1/do.call(dfun, c(list(x = qpop), params))
    }

  } else stop("'method' must be either 'qor' or 'density'.\n")

  return(out)

}

