#' qor.ln
#'
#' calculates the QOR (Quantile optimality ratio) for the log normal
#'
#' @param u a numeric vector of probability values in \[0, 1\] indicating all quantiles to be estimated.
#' @return QOR function for the log-normal
#' @noRd

qor.ln <- function(u){
  # QOR function for the log-normal
  q.n.0 <- 1/dnorm(qnorm(u))
  q.n.1 <- qnorm(u)*q.n.0^2
  q.n.2 <- (1 + 2*qnorm(u)^2)*q.n.0^3
  1/(q.n.0^2 + 3*q.n.1 + q.n.2/q.n.0)
}
