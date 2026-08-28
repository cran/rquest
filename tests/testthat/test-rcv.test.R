test_that("rcv.test returns htest class for MAD and IQR versions", {

  set.seed(1234)
  x <- rnorm(100, mean = 8)

  expect_s3_class(
    rcv.test(x, numerator = "mad"),
    "htest"
  )

  expect_s3_class(
    rcv.test(x, numerator = "iqr"),
    "htest"
  )
})


test_that("rcv.test works for two samples", {

  set.seed(1234)
  x <- rnorm(100, mean = 8)
  y <- rnorm(120, mean = 10)

  expect_s3_class(
    rcv.test(x, y, numerator = "mad"),
    "htest"
  )

  expect_s3_class(
    rcv.test(x, y, numerator = "iqr"),
    "htest"
  )
})


test_that("rcv.test removes NA values with warnings", {

  set.seed(1234)

  x <- c(rnorm(100, mean = 8), NA)
  y <- rnorm(120, mean = 10)

  expect_warning(
    rcv.test(x, numerator = "mad"),
    "missing values removed"
  )

  expect_warning(
    rcv.test(y, x, numerator = "mad"),
    "missing values removed"
  )
})


test_that("rcv.test validates x and y", {

  expect_error(
    rcv.test("abc"),
    "Argument 'x' must be numeric.",
    fixed = TRUE
  )

  expect_error(
    rcv.test(1:10, y = "abc"),
    "argument 'y' must be numeric",
    fixed = TRUE
  )
})


test_that("rcv.test validates numerator", {

  expect_error(
    rcv.test(
      1:20,
      numerator = "wrong"
    ),
    "Argument 'numerator' must be either 'mad' or 'iqr'.",
    fixed = TRUE
  )
})


test_that("MAD rCV estimate matches manual calculation", {

  set.seed(1234)
  x <- rnorm(100, mean = 8)

  res <- rcv.test(
    x,
    numerator = "mad",
    log.transf = FALSE,
    back.transf = FALSE
  )

  expected <- 1.4826 *
    mad(x, constant = 1) /
    median(x)

  expect_equal(
    as.numeric(res$estimate),
    as.numeric(expected)
  )
})


test_that("IQR rCV estimate matches q.test", {

  set.seed(1234)
  x <- rnorm(100, mean = 8)

  res.rcv <- rcv.test(
    x,
    numerator = "iqr"
  )

  res.q <- q.test(
    x,
    measure = "rCViqr",
    log.transf = TRUE,
    back.transf = TRUE
  )

  expect_equal(
    as.numeric(res.rcv$estimate),
    as.numeric(res.q$estimate)
  )

  expect_equal(
    res.rcv$conf.int,
    res.q$conf.int
  )
})


test_that("rcv.test passes alternative to q.test for IQR version", {

  set.seed(1234)
  x <- rnorm(100, mean = 8)

  res.less <- rcv.test(
    x,
    numerator = "iqr",
    alternative = "less"
  )

  res.greater <- rcv.test(
    x,
    numerator = "iqr",
    alternative = "greater"
  )

  expect_equal(
    res.less$alternative,
    "less"
  )

  expect_equal(
    res.greater$alternative,
    "greater"
  )
})


test_that("MAD less alternative gives valid one-sided confidence interval", {

  set.seed(1234)
  x <- rnorm(100, mean = 8)

  res <- rcv.test(
    x,
    numerator = "mad",
    alternative = "less"
  )

  expect_equal(
    unname(res$conf.int[1]),
    0
  )

  expect_true(
    res$conf.int[1] <= res$conf.int[2]
  )
})


test_that("MAD greater alternative gives valid one-sided confidence interval", {

  set.seed(1234)
  x <- rnorm(100, mean = 8)

  res <- rcv.test(
    x,
    numerator = "mad",
    alternative = "greater"
  )

  expect_true(
    is.infinite(res$conf.int[2])
  )

  expect_true(
    res$conf.int[1] <= res$conf.int[2]
  )
})


test_that("MAD two-sided confidence interval contains estimate", {

  set.seed(1234)
  x <- rnorm(100, mean = 8)

  res <- rcv.test(
    x,
    numerator = "mad",
    alternative = "two.sided"
  )

  expect_true(
    res$conf.int[1] <= res$estimate
  )

  expect_true(
    res$estimate <= res$conf.int[2]
  )
})


test_that("two-sample MAD back-transformed estimate is ratio of rCVs", {

  set.seed(1234)
  x <- rnorm(100, mean = 8)
  y <- rnorm(120, mean = 10)

  res <- rcv.test(
    x,
    y,
    numerator = "mad",
    log.transf = TRUE,
    back.transf = TRUE
  )

  rcv.x <- 1.4826 *
    mad(x, constant = 1) /
    median(x)

  rcv.y <- 1.4826 *
    mad(y, constant = 1) /
    median(y)

  expected <- rcv.x / rcv.y

  expect_equal(
    as.numeric(res$estimate),
    as.numeric(expected)
  )
})


test_that("two-sample MAD without log transformation gives difference in rCVs", {

  set.seed(1234)
  x <- rnorm(100, mean = 8)
  y <- rnorm(120, mean = 10)

  res <- rcv.test(
    x,
    y,
    numerator = "mad",
    log.transf = FALSE,
    back.transf = FALSE
  )

  rcv.x <- 1.4826 *
    mad(x, constant = 1) /
    median(x)

  rcv.y <- 1.4826 *
    mad(y, constant = 1) /
    median(y)

  expected <- rcv.x - rcv.y

  expect_equal(
    as.numeric(res$estimate),
    as.numeric(expected)
  )
})


test_that("IQR version supports both variance methods", {

  set.seed(1234)
  x <- rlnorm(100)

  res.qor <- rcv.test(
    x,
    numerator = "iqr",
    var.method = "qor"
  )

  res.density <- rcv.test(
    x,
    numerator = "iqr",
    var.method = "density"
  )

  expect_equal(
    as.numeric(res.qor$estimate),
    as.numeric(res.density$estimate)
  )

  expect_false(
    isTRUE(
      all.equal(
        res.qor$conf.int,
        res.density$conf.int
      )
    )
  )
})


test_that("rcv.test passes density arguments in MAD branch", {

  set.seed(1234)
  x <- rnorm(100, mean = 8)

  expect_s3_class(
    rcv.test(
      x,
      numerator = "mad",
      bw = 0.5
    ),
    "htest"
  )
})


test_that("rcv.test passes density arguments to both samples", {

  set.seed(1234)
  x <- rnorm(100, mean = 8)
  y <- rnorm(120, mean = 10)

  expect_s3_class(
    rcv.test(
      x,
      y,
      numerator = "mad",
      bw = 0.5
    ),
    "htest"
  )
})


test_that("log transformation requires positive estimates", {

  x <- c(-10, -9, -8, -7, -6, -5, -4, -3, -2, -1)

  expect_error(
    rcv.test(
      x,
      numerator = "mad",
      log.transf = TRUE
    ),
    "Estimates must be positive to use the log transformation.",
    fixed = TRUE
  )
})


test_that("rcv.test output labels are correct", {

  set.seed(1234)
  x <- rnorm(100, mean = 8)
  y <- rnorm(120, mean = 10)

  one <- rcv.test(
    x,
    numerator = "mad"
  )

  two <- rcv.test(
    x,
    y,
    numerator = "mad"
  )

  expect_equal(
    one$method,
    "One sample test of the robust coefficient of variation (MAD/median)"
  )

  expect_equal(
    two$method,
    "Two sample test of the robust coefficient of variation (MAD/median)"
  )

  expect_match(
    names(one$estimate),
    "Robust CV"
  )

  expect_match(
    names(two$estimate),
    "ratio of Robust CVs"
  )
})
