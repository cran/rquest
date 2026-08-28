test_that("qcov returns a matrix", {
  expect_true(
    is.matrix(qcov(1:10, c(0.25, 0.5, 0.75)))
  )
})


test_that("qcov gives consistent variances for repeated quantiles", {

  V <- qcov(1:10, c(0.25, 0.25))
  Vsingle <- qcov(1:10, 0.25)

  expect_equal(V[1, 1], V[1, 2])
  expect_equal(V[1, 1], V[2, 1])
  expect_equal(V[1, 1], V[2, 2])
  expect_equal(V[1, 1], Vsingle[1, 1])
})


test_that("qcov works with density method", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  V <- qcov(
    x,
    u,
    method = "density"
  )

  expect_true(is.matrix(V))
  expect_equal(dim(V), c(3, 3))
  expect_true(isSymmetric(V))
  expect_true(all(diag(V) >= 0))
})


test_that("qcov density method works with supplied distribution parameters", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  V <- qcov(
    x,
    u,
    method = "density",
    dist = "norm",
    params = list(mean = 0, sd = 1)
  )

  expect_true(is.matrix(V))
  expect_equal(dim(V), c(3, 3))
  expect_true(isSymmetric(V))
  expect_true(all(is.finite(V)))
})


test_that("qcov with known normal parameters matches manual calculation", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  qd <- 1 / dnorm(
    qnorm(u, mean = 0, sd = 1),
    mean = 0,
    sd = 1
  )

  u1u <- u %*% t(1 - u)
  u1u <- pmin(u1u, t(u1u))

  expected <- u1u * tcrossprod(qd) / length(x)

  observed <- qcov(
    x,
    u,
    method = "density",
    dist = "norm",
    params = list(mean = 0, sd = 1)
  )

  expect_equal(
    unname(observed),
    unname(expected)
  )
})


test_that("qcov validates inputs", {

  expect_error(
    qcov("abc", c(0.25, 0.5)),
    "must be numeric"
  )

  expect_error(
    qcov(1:10, c(0, 0.5)),
    "between, but not including, 0 and 1"
  )

  expect_error(
    qcov(1:10, c(0.25, NA)),
    "between, but not including, 0 and 1"
  )

  expect_error(
    qcov(1:10, c(0.25, 0.5), method = "wrong"),
    "must be either 'qor' or 'density'"
  )
})
