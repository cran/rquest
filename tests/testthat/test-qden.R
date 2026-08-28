test_that("qden returns numeric output of correct length", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  out <- qden(x, u, method = "density")

  expect_true(is.numeric(out))
  expect_length(out, length(u))
})


test_that("qden density method works", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  out <- qden(
    x,
    u,
    method = "density"
  )

  expect_true(all(is.finite(out)))
  expect_true(all(out > 0))
})


test_that("qden density method with supplied normal parameters matches manual calculation", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  observed <- qden(
    x,
    u,
    dist = "norm",
    method = "density",
    params = list(mean = 0, sd = 1)
  )

  qpop <- qnorm(
    u,
    mean = 0,
    sd = 1
  )

  expected <- 1 / dnorm(
    qpop,
    mean = 0,
    sd = 1
  )

  expect_equal(
    observed,
    expected
  )
})


test_that("qden qor method works with normal distribution", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  out <- qden(
    x,
    u,
    dist = "norm",
    method = "qor"
  )

  expect_true(is.numeric(out))
  expect_length(out, length(u))
  expect_true(all(is.finite(out)))
})


test_that("qden default GLD qor method works", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  out <- qden(
    x,
    u,
    method = "qor"
  )

  expect_true(is.numeric(out))
  expect_length(out, length(u))
  expect_true(all(is.finite(out)))
})


test_that("qden supports all kernels", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  kernels <- c(
    "epanechnikov",
    "gaussian",
    "rectangular",
    "triangular",
    "biweight",
    "cosine",
    "optcosine"
  )

  for (k in kernels) {

    out <- qden(
      x,
      u,
      dist = "norm",
      method = "qor",
      kernel = k
    )

    expect_true(is.numeric(out))
    expect_length(out, length(u))
    expect_true(all(is.finite(out)))
  }
})


test_that("qden rejects invalid kernel", {

  set.seed(1234)
  x <- rnorm(100)

  expect_error(
    qden(
      x,
      c(0.25, 0.5),
      dist = "norm",
      method = "qor",
      kernel = "wrong"
    ),
    "should be one of"
  )
})


test_that("qden validates method", {

  expect_error(
    qden(
      1:10,
      c(0.25, 0.5),
      method = "wrong"
    ),
    "'method' must be either 'qor' or 'density'.",
    fixed = TRUE
  )
})


test_that("qden validates x", {

  expect_error(
    qden(
      "abc",
      c(0.25, 0.5)
    ),
    "Argument 'x' must be numeric.",
    fixed = TRUE
  )
})


test_that("qden validates u", {

  expect_error(
    qden(
      1:10,
      c(0, 0.5),
      method = "density"
    ),
    "between, but not including, 0 and 1"
  )

  expect_error(
    qden(
      1:10,
      c(0.25, 1),
      method = "density"
    ),
    "between, but not including, 0 and 1"
  )

  expect_error(
    qden(
      1:10,
      c(0.25, NA),
      method = "density"
    ),
    "between, but not including, 0 and 1"
  )
})


test_that("qden large sample switch gives warning", {

  set.seed(1234)
  x <- rlnorm(1001)

  expect_warning(
    qden(
      x,
      c(0.25, 0.5, 0.75),
      method = "qor"
    ),
    "lnorm"
  )
})


test_that("qden large sample switch can be disabled", {

  set.seed(1234)
  x <- rlnorm(1001)

  expect_no_warning(
    qden(
      x,
      c(0.25, 0.5, 0.75),
      method = "qor",
      large.n.switch = FALSE
    )
  )
})
