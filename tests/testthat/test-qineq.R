test_that("qineq returns htest class", {

  expect_s3_class(
    qineq(1:10),
    "htest"
  )
})


test_that("qineq removes NA values with warnings", {

  x <- c(1:20, NA)
  y <- 2:21

  expect_warning(
    qineq(x),
    "missing values removed"
  )

  expect_warning(
    qineq(y, x),
    "missing values removed"
  )
})


test_that("qineq validates measure argument", {

  expect_error(
    qineq(1:10, measure = "Gini"),
    "Unknown inequality measure specified"
  )
})


test_that("qineq output is correct for standard cases", {

  expect_equal(
    qineq(1:10)$method,
    "One sample test of the QRI statistic"
  )

  expect_equal(
    qineq(1:10, 1:10)$method,
    "Two sample test of the QRI statistic"
  )

  expect_equal(
    qineq(1:10, measure = "G2")$method,
    "One sample test of the G2 statistic"
  )

  expect_equal(
    qineq(1:10)$data.name,
    "1:10"
  )

  expect_equal(
    as.numeric(qineq(1:10, 1:10)$statistic),
    0
  )

  expect_equal(
    qineq(1:10, alternative = "less")$alternative,
    "less"
  )

  expect_equal(
    as.numeric(qineq(1:10, 1:10)$estimate),
    0
  )

  expect_equal(
    as.numeric(qineq(1:10, 1:10, true.ineq = 5)$null.value),
    5
  )
})


test_that("all built-in inequality measures work", {

  x <- 1:20

  measures <- c("QRI", "G1", "G2", "G3", "S1", "S2")

  for (m in measures) {
    expect_s3_class(
      qineq(x, measure = m),
      "htest"
    )
  }
})


test_that("formula-defined inequality measure works", {

  x <- 1:20

  res <- qineq(
    x,
    measure = ~ Q(p / 2) / Q(1 - p / 2)
  )

  expect_s3_class(res, "htest")
  expect_match(res$method, "user defined statistic")
})


test_that("list-defined inequality measure works", {

  x <- 1:20
  J <- 20
  p <- (1:J - 0.5) / J

  u <- sort(c(p / 2, 1 - p / 2))

  num <- cbind(
    diag(rep(1, J)),
    matrix(0, J, J)
  )

  den <- cbind(
    matrix(0, J, J),
    diag(J)[, J:1]
  )

  res <- qineq(
    x,
    J = J,
    measure = list(
      u = u,
      coef1 = num,
      coef2 = den
    )
  )

  expect_s3_class(res, "htest")
  expect_match(res$method, "user defined statistic")
})


test_that("qineq supports both variance methods", {

  set.seed(1234)
  x <- rlnorm(100)

  res.qor <- qineq(
    x,
    var.method = "qor"
  )

  res.density <- qineq(
    x,
    var.method = "density"
  )

  expect_equal(
    as.numeric(res.qor$estimate),
    as.numeric(res.density$estimate)
  )

  expect_false(
    isTRUE(all.equal(
      as.numeric(res.qor$statistic),
      as.numeric(res.density$statistic)
    ))
  )
})


test_that("two-sample less confidence interval uses valid lower bound", {

  set.seed(1234)

  x <- rlnorm(500, meanlog = 0, sdlog = 0.2)
  y <- rlnorm(500, meanlog = 0, sdlog = 1)

  res <- qineq(
    x,
    y,
    measure = "QRI",
    alternative = "less"
  )

  expect_equal(
    unname(res$conf.int[1]),
    -1
  )

  expect_true(
    res$conf.int[1] <= res$conf.int[2]
  )
})
