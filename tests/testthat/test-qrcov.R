test_that("qrcov returns expected list structure", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  coef1 <- matrix(
    c(0, 0, 1,
      1, 0, 0),
    nrow = 2,
    byrow = TRUE
  )

  coef2 <- matrix(
    c(0, 1, 0,
      0, 1, 0),
    nrow = 2,
    byrow = TRUE
  )

  res <- qrcov(
    x,
    u,
    coef1 = coef1,
    coef2 = coef2
  )

  expect_type(res, "list")
  expect_named(res, c("ratios", "cov"))
  expect_length(res$ratios, 2)
  expect_true(is.matrix(res$cov))
  expect_equal(dim(res$cov), c(2, 2))
})


test_that("qrcov ratio estimates match manual calculations", {

  x <- 1:20
  u <- c(0.25, 0.5, 0.75)

  coef1 <- matrix(
    c(0, 0, 1,
      1, 0, 0),
    nrow = 2,
    byrow = TRUE
  )

  coef2 <- matrix(
    c(0, 1, 0,
      0, 1, 0),
    nrow = 2,
    byrow = TRUE
  )

  res <- qrcov(
    x,
    u,
    coef1 = coef1,
    coef2 = coef2
  )

  qest <- quantile(
    x,
    u,
    type = 8
  )

  expected <- c(
    qest[3] / qest[2],
    qest[1] / qest[2]
  )

  expect_equal(
    unname(res$ratios),
    unname(expected)
  )
})


test_that("qrcov covariance matrix is symmetric", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  coef1 <- matrix(
    c(0, 0, 1,
      1, 0, 0),
    nrow = 2,
    byrow = TRUE
  )

  coef2 <- matrix(
    c(0, 1, 0,
      0, 1, 0),
    nrow = 2,
    byrow = TRUE
  )

  res <- qrcov(
    x,
    u,
    coef1 = coef1,
    coef2 = coef2
  )

  expect_true(
    isSymmetric(res$cov)
  )

  expect_true(
    all(diag(res$cov) >= 0)
  )
})


test_that("qrcov covariance matrix has correct names", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  coef1 <- matrix(
    c(0, 0, 1,
      1, 0, 0),
    nrow = 2,
    byrow = TRUE
  )

  coef2 <- matrix(
    c(0, 1, 0,
      0, 1, 0),
    nrow = 2,
    byrow = TRUE
  )

  res <- qrcov(
    x,
    u,
    coef1 = coef1,
    coef2 = coef2
  )

  expect_equal(
    rownames(res$cov),
    c("R1", "R2")
  )

  expect_equal(
    colnames(res$cov),
    c("R1", "R2")
  )
})


test_that("qrcov works with density method", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  coef1 <- matrix(
    c(0, 0, 1,
      1, 0, 0),
    nrow = 2,
    byrow = TRUE
  )

  coef2 <- matrix(
    c(0, 1, 0,
      0, 1, 0),
    nrow = 2,
    byrow = TRUE
  )

  res <- qrcov(
    x,
    u,
    coef1 = coef1,
    coef2 = coef2,
    method = "density"
  )

  expect_true(is.numeric(res$ratios))
  expect_true(is.matrix(res$cov))
  expect_true(all(is.finite(res$cov)))
})


test_that("qrcov works with qor method", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  coef1 <- matrix(
    c(0, 0, 1,
      1, 0, 0),
    nrow = 2,
    byrow = TRUE
  )

  coef2 <- matrix(
    c(0, 1, 0,
      0, 1, 0),
    nrow = 2,
    byrow = TRUE
  )

  res <- qrcov(
    x,
    u,
    coef1 = coef1,
    coef2 = coef2,
    method = "qor"
  )

  expect_true(is.numeric(res$ratios))
  expect_true(is.matrix(res$cov))
  expect_true(all(is.finite(res$cov)))
})


test_that("qrcov gives same ratios for qor and density methods", {

  set.seed(1234)
  x <- rnorm(100)
  u <- c(0.25, 0.5, 0.75)

  coef1 <- matrix(
    c(0, 0, 1,
      1, 0, 0),
    nrow = 2,
    byrow = TRUE
  )

  coef2 <- matrix(
    c(0, 1, 0,
      0, 1, 0),
    nrow = 2,
    byrow = TRUE
  )

  res.qor <- qrcov(
    x,
    u,
    coef1 = coef1,
    coef2 = coef2,
    method = "qor"
  )

  res.density <- qrcov(
    x,
    u,
    coef1 = coef1,
    coef2 = coef2,
    method = "density"
  )

  expect_equal(
    res.qor$ratios,
    res.density$ratios
  )
})


test_that("qrcov validates x", {

  u <- c(0.25, 0.5, 0.75)

  coef1 <- matrix(
    c(0, 0, 1),
    nrow = 1
  )

  coef2 <- matrix(
    c(0, 1, 0),
    nrow = 1
  )

  expect_error(
    qrcov(
      "abc",
      u,
      coef1 = coef1,
      coef2 = coef2
    ),
    "Argument 'x' must be numeric.",
    fixed = TRUE
  )
})


test_that("qrcov validates u", {

  coef1 <- matrix(
    c(0, 0, 1),
    nrow = 1
  )

  coef2 <- matrix(
    c(0, 1, 0),
    nrow = 1
  )

  expect_error(
    qrcov(
      1:20,
      c(0, 0.5, 0.75),
      coef1 = coef1,
      coef2 = coef2
    ),
    "between, but not including, 0 and 1"
  )

  expect_error(
    qrcov(
      1:20,
      c(0.25, NA, 0.75),
      coef1 = coef1,
      coef2 = coef2
    ),
    "between, but not including, 0 and 1"
  )
})


test_that("qrcov requires numeric matrix coefficients", {

  u <- c(0.25, 0.5, 0.75)

  expect_error(
    qrcov(
      1:20,
      u,
      coef1 = c(0, 0, 1),
      coef2 = matrix(c(0, 1, 0), nrow = 1)
    ),
    "Arguments 'coef1' and 'coef2' both must be a numeric matrix.",
    fixed = TRUE
  )

  expect_error(
    qrcov(
      1:20,
      u,
      coef1 = matrix(c("a", "b", "c"), nrow = 1),
      coef2 = matrix(c(0, 1, 0), nrow = 1)
    ),
    "Arguments 'coef1' and 'coef2' both must be a numeric matrix.",
    fixed = TRUE
  )
})


test_that("qrcov requires matching coefficient dimensions", {

  u <- c(0.25, 0.5, 0.75)

  coef1 <- matrix(
    c(0, 0, 1,
      1, 0, 0),
    nrow = 2,
    byrow = TRUE
  )

  coef2 <- matrix(
    c(0, 1, 0),
    nrow = 1
  )

  expect_error(
    qrcov(
      1:20,
      u,
      coef1 = coef1,
      coef2 = coef2
    ),
    "same dimension"
  )
})


test_that("qrcov requires number of columns to equal length of u", {

  u <- c(0.25, 0.5, 0.75)

  coef1 <- matrix(
    c(1, 0),
    nrow = 1
  )

  coef2 <- matrix(
    c(0, 1),
    nrow = 1
  )

  expect_error(
    qrcov(
      1:20,
      u,
      coef1 = coef1,
      coef2 = coef2
    ),
    "number of columns equal to the length of argument 'u'"
  )
})
