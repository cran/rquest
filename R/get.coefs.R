#' Extract coefficients for a quantile-ratio formula
#'
#' Internal helper function used to convert a user-supplied quantile-ratio
#' formula into the quantile probabilities and coefficient matrices required
#' by `qineq()`.
#'
#' @param formula A formula defining a ratio of linear combinations of
#'   quantiles.
#' @param pvec A numeric vector of probability values.
#'
#' @return A list containing `u`, `coef1`, and `coef2`.
#'
#' @keywords internal
#' @noRd

get.coefs <- function(formula, pvec) {
  # length of formula should be 2.  With formula[[1]] = `~`, formula[[2]] = equation

  eq <- formula[[2]] # gets the form of the formula
  if(eq[[1]] != as.name("/")) stop("Please make sure your formula is a valid ratio expression using `/`.\n")

  # length of eq is 3.  eq[[1]] is the middle operator, should be `/` tested above
  # eq[[2]] is the numerator (LHS) of `/` and eq[[3]] the denominator
  num <- eq[[2]]
  den <- eq[[3]]

  get.u <- function(x, p) {
    if(is.call(x)){
      if(x[[1]] == as.name("Q")) return(as.numeric(eval(x[[2]], list(p = p))))
      else  return(unlist(lapply(as.list(x)[-1], get.u, p = p), use.names = FALSE))
    } else numeric(0)
  }

  all.us <- unlist(lapply(pvec,function(p) c(get.u(num, p), get.u(den, p))))
  uvec <- sort(unique(all.us))
  len.u <- length(uvec)

  coef1 <- matrix(0, length(pvec), len.u)
  coef2 <- coef1

  lin <- function(x, p) {
    if (is.numeric(x)){
      # finds a constant scalar multiplier: e.g. 2*
      return(list(cst = as.numeric(x), q = rep(0, len.u)))
    } else if (is.symbol(x)) {
      # for finding the p* multipliers
      if (as.character(x) != "p")
        stop("The only non-constant variable allowed outside of `Q(...)` is 'p'.\n",
             call. = FALSE)
      return(list(cst = p, q = rep(0, len.u)))
    } else if (!is.call(x)) stop("There is something wrong with your formula.  Please check the help files.\n",
                                 call. = FALSE)

    # Depends on x.
    # If x = Q(p/2), x[[1]] = `Q`, x[[2]] = p/2
    # If x = Q(p/2) + Q(1-p/2), x[[1]] = `+`, x[[2]] = Q(p/2) and x[[3]] = Q(1 - p/2)
    # If x = 2*p*Q(p/2), x[[1]] = `*`, x[[2]] = 2*p, X[[3]] = Q(p/2)
    op <- as.character(x[[1]]) # should be 'Q', '(', '*', '+', '-'
    if (op == "Q") {
      u <- as.numeric(eval(x[[2]], list(p = p)))
      j <- match(round(u, 12), round(uvec, 12))
      q <- rep(0, len.u); q[j] <- 1
      return(list(cst = 0, q = q))
    }

    a <- lin(x[[2]], p) # LHS of operator
    if(length(x) > 2) b <- lin(x[[3]], p) # RHS of operator
    else b <- NULL

    if (op == "+") return(list(cst = a$cst + b$cst, q = a$q + b$q))
    if (op == "-") return(list(cst = a$cst - b$cst, q = a$q - b$q))
    # below gets LHS of * and RHS of *
    if (op == "*") {
      if (any(a$q != 0) && any(b$q != 0)) stop("Formula cannot contain a multiplication of two Q(..)'s.\n",
                                               call. = FALSE)
      if (any(a$q != 0)) return(list(cst = 0, q = a$q * b$cst))
      if (any(b$q != 0)) return(list(cst = 0, q = b$q * a$cst))
      return(list(cst = a$cst * b$cst, q = rep(0, len.u)))
    }
    if (op == "(") return(a)

    stop("Operator `", as.character(op), "` not allowed. Please check the help files.\n",
         call. = FALSE)
  }

  for (i in seq_along(pvec)) {
    coef1[i, ] <- lin(num, pvec[i])$q
    coef2[i, ] <- lin(den, pvec[i])$q
  }

  list(u = uvec, coef1 = coef1, coef2 = coef2)
}
