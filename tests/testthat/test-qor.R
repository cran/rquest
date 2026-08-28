test_that("qor returns a list with expected components", {

  u <- c(0.25, 0.5, 0.75)

  res <- qor(
    u,
    dist = "norm"
  )

  expect_type(res, "list")
  expect_named(res, c("qor", "params"))
  expect_length(res$qor, length(u))
})


test_that("qor normal distribution matches closed-form expression", {

  u <- c(0.25, 0.5, 0.75)

  observed <- qor(
    u,
    dist = "norm"
  )$qor

  zu <- qnorm(u)

  expected <- dnorm(zu)^2 / (1 + 2 * zu^2)

  expect_equal(
    observed,
    expected
  )
})


test_that("qor exponential distribution matches closed-form expression", {

  u <- c(0.25, 0.5, 0.75)

  observed <- qor(
    u,
    dist = "exp"
  )$qor

  expected <- (1 - u)^2 / 2

  expect_equal(
    observed,
    expected
  )
})


test_that("qor lognormal distribution works with x", {

  set.seed(1234)
  x <- rlnorm(100)
  u <- c(0.25, 0.5, 0.75)

  res <- qor(
    u,
    dist = "lnorm",
    x = x
  )

  expect_true(is.numeric(res$qor))
  expect_length(res$qor, length(u))
  expect_true(all(is.finite(res$qor)))
})


test_that("qor lognormal distribution matches manual calculation with supplied parameters", {

  u <- c(0.25, 0.5, 0.75)
  sdlog <- 0.8

  observed <- qor(
    u,
    dist = "lnorm",
    params = list(sdlog = sdlog)
  )$qor

  zu <- qnorm(u)

  expected <- dnorm(zu)^2 /
    (1 + sdlog^2 + 3 * sdlog * zu + 2 * zu^2)

  expect_equal(
    observed,
    expected
  )
})


test_that("qor GLD works with fitted parameters from x", {

  set.seed(1234)
  x <- rlnorm(100)
  u <- c(0.25, 0.5, 0.75)

  res <- qor(
    u,
    dist = "gl",
    x = x
  )

  expect_true(is.numeric(res$qor))
  expect_length(res$qor, length(u))
  expect_true(all(is.finite(res$qor)))
})


test_that("qor GLD works with supplied parameters", {

  u <- c(0.25, 0.5, 0.75)

  params <- list(
    lambda1 = 0,
    lambda2 = 1,
    lambda3 = 0.2,
    lambda4 = 0.2
  )

  res <- qor(
    u,
    dist = "gl",
    params = params
  )

  expect_true(is.numeric(res$qor))
  expect_length(res$qor, length(u))
  expect_true(all(is.finite(res$qor)))
})


test_that("qor user-supplied quantile function works", {

  u <- c(0.2, 0.5, 0.8)

  res <- qor(
    u,
    dist = qbeta,
    params = list(
      shape1 = 2,
      shape2 = 5
    )
  )

  expect_true(is.numeric(res$qor))
  expect_length(res$qor, length(u))
  expect_true(all(is.finite(res$qor)))
})


test_that("qor user-supplied quantile function warns when params are NULL", {

  u <- c(0.25, 0.5, 0.75)

  expect_warning(
    qor(
      u,
      dist = qnorm
    ),
    "Using the default parameter values"
  )
})


test_that("qor requires x or params for GLD", {

  expect_error(
    qor(
      c(0.25, 0.5, 0.75),
      dist = "gl"
    ),
    "Must provide 'x' if 'params' is NULL.",
    fixed = TRUE
  )
})


test_that("qor requires x or params for lognormal distribution", {

  expect_error(
    qor(
      c(0.25, 0.5, 0.75),
      dist = "lnorm"
    ),
    "Must provide 'x' if 'params' is NULL.",
    fixed = TRUE
  )
})


test_that("qor rejects non-positive x for lognormal distribution", {

  expect_error(
    qor(
      c(0.25, 0.5, 0.75),
      dist = "lnorm",
      x = c(-1, 1, 2, 3)
    ),
    "Values in 'x' must be positive if using 'lnorm'.",
    fixed = TRUE
  )
})


test_that("qor validates u", {

  expect_error(
    qor(
      c(0, 0.5),
      dist = "norm"
    ),
    "between, but not including, 0 and 1"
  )

  expect_error(
    qor(
      c(0.25, 1),
      dist = "norm"
    ),
    "between, but not including, 0 and 1"
  )

  expect_error(
    qor(
      c(0.25, NA),
      dist = "norm"
    ),
    "between, but not including, 0 and 1"
  )
})


test_that("qor rejects unknown named distribution", {

  expect_error(
    qor(
      c(0.25, 0.5),
      dist = "wrong"
    ),
    "Unknown argument for 'dist'.",
    fixed = TRUE
  )
})


test_that("qor rejects invalid dist type", {

  expect_error(
    qor(
      c(0.25, 0.5),
      dist = 123
    ),
    "'dist' must be a character string or a function.",
    fixed = TRUE
  )
})
