test_that("q.test returns htest class", {
  expect_s3_class(q.test(c(1:10)), "htest")
})


test_that("q.test removes NA values with warnings", {

  x <- c(1:20, NA)
  y <- 2:21

  expect_warning(
    q.test(x),
    "missing values removed"
  )

  expect_warning(
    q.test(y, x),
    "missing values removed"
  )
})


test_that("ratio warning and log transformation behaviour", {

  expect_warning(
    q.test(
      c(1:10),
      u = c(0.25, 0.5),
      coef = c(-1, 1),
      u2 = c(0.5, 0.75),
      coef2 = c(-1, 1)
    ),
    "consider using a log transformation"
  )

  expect_no_warning(
    q.test(
      c(1:10),
      u = c(0.25, 0.5),
      coef = c(-1, 1),
      u2 = c(0.5, 0.75),
      coef2 = c(-1, 1),
      log.transf = TRUE
    )
  )
})


test_that("q.test validates arguments", {

  expect_error(
    q.test(c(1:10), coef2 = 1),
    "When using u2 and coef2, you also need to specify both u and coef."
  )

  expect_error(
    q.test(c(1:10), coef = 1, coef2 = 1),
    "When using u2 and coef2, you also need to specify both u and coef."
  )

  expect_error(
    q.test(c(1:10), measure = "bowley", p = 0.5),
    "Argument p must be a numeric value in (0, 1) except 0.5.",
    fixed = TRUE
  )

  expect_error(
    q.test(c(1:10), measure = "measure"),
    "Unknown choice for measure."
  )

  expect_error(
    q.test(c(1:10), measure = "groen"),
    "Use one of 'groenR' or 'groenL'."
  )

  expect_error(
    q.test(c(1:10), measure = "groenL", p = 0.5),
    "Argument p must be a numeric value in (0, 1) except 0.5.",
    fixed = TRUE
  )

  expect_error(
    q.test(c(1:10), measure = "lqw", p = 0.75),
    "Argument p must be a numeric value in (0, 1/2).",
    fixed = TRUE
  )

  expect_error(
    q.test(c(1:10), measure = "rqw", p = 0.25),
    "Argument p must be a numeric value in (1/2, 1).",
    fixed = TRUE
  )

  expect_error(
    q.test(c(1:10), measure = "qr0.250.75"),
    "For quantile ratios, measure must be in format qrxxyy where xx and yy are integer numbers."
  )

  expect_error(
    q.test(
      c(1:10),
      u = c(0.25, 0.75),
      coef = c(1, "a")
    ),
    "Probability and coefficient arguments (i.e. 'u', 'coef', 'u2', 'coef2') must be either numeric or not specified (default to NULL).",
    fixed = TRUE
  )

  expect_error(
    q.test(c(1:10), coef = 1),
    "Argument u required if coef is a numeric vector."
  )

  expect_error(
    q.test(c(1:10), u = 0.25, coef = c(1, -1)),
    "Length of u needs to be equal to the length of coef."
  )

  expect_error(
    q.test(
      c(1:10),
      u = c(0.25, 0.75),
      coef = matrix(c(1, "a", 0, 1), 2, 2)
    ),
    "Probability and coefficient arguments (i.e. 'u', 'coef', 'u2', 'coef2') must be either numeric or not specified (default to NULL).",
    fixed = TRUE
  )

  expect_error(
    q.test(
      c(1:10),
      coef = matrix(c(1, 0, 0, 1), 2, 2)
    ),
    "Argument u required if coef is a numeric matrix."
  )

  expect_error(
    q.test(
      c(1:10),
      u = c(0.25, 0.75),
      coef = matrix(c(1, 0, 1, 0, 1, 1), 3, 2)
    ),
    "Matrix coef needs to have dimensions: ncol(coef)=length(u) and nrow(coef)=2.",
    fixed = TRUE
  )

  expect_error(
    q.test(c(1:10), var.method = "method"),
    "Argument var.method must be either 'qor' or 'density'.",
    fixed = TRUE
  )
})


test_that("LCQ calculation is correct", {

  expect_equal(
    as.numeric(
      suppressWarnings(
        q.test(
          c(1:10),
          u = c(0.25, 0.75),
          coef = c(-2, 4),
          u2 = c(0.5, 0.75),
          coef2 = c(-2, 3)
        )
      )$estimate
    ),
    2
  )
})


test_that("vector and matrix ratio inputs are equivalent", {

  res1 <- suppressWarnings(
    q.test(
      c(1:10),
      u = c(0.25, 0.75),
      coef = c(-1, 1),
      u2 = c(0.5, 0.75),
      coef2 = c(-1, 1)
    )
  )

  res2 <- suppressWarnings(
    q.test(
      c(1:10),
      u = c(0.25, 0.5, 0.75),
      coef = matrix(
        c(-1, 0, 1,
          0, -1, 1),
        2, 3,
        byrow = TRUE
      )
    )
  )

  expect_equal(res1, res2)
})


test_that("rCViqr predefined and user-defined versions agree", {

  res1 <- suppressWarnings(
    q.test(c(1:10), measure = "rCViqr")
  )

  res2 <- suppressWarnings(
    q.test(
      c(1:10),
      u = c(0.25, 0.75),
      coef = c(-0.75, 0.75),
      u2 = 0.5
    )
  )

  expect_equal(res1$conf.int, res2$conf.int)
})


test_that("IQR predefined and user-defined versions agree", {

  expect_equal(
    q.test(c(1:10), measure = "IQR")$conf.int,
    q.test(
      c(1:10),
      u = c(0.25, 0.75),
      coef = c(-1, 1)
    )$conf.int
  )
})


test_that("median predefined and user-defined versions agree", {

  expect_equal(
    q.test(c(1:10), measure = "median")$conf.int,
    q.test(c(1:10), u = 0.5)$conf.int
  )
})


test_that("Bowley predefined and user-defined versions agree", {

  res1 <- suppressWarnings(
    q.test(c(1:10), measure = "bowley")
  )

  res2 <- suppressWarnings(
    q.test(
      c(1:10),
      u = c(0.25, 0.5, 0.75),
      coef = c(1, -2, 1),
      u2 = c(0.25, 0.75),
      coef2 = c(-1, 1)
    )
  )

  expect_equal(res1$conf.int, res2$conf.int)
})


test_that("Kelly predefined and user-defined versions agree", {

  res1 <- suppressWarnings(
    q.test(c(1:10), measure = "kelly")
  )

  res2 <- suppressWarnings(
    q.test(
      c(1:10),
      u = c(0.1, 0.5, 0.9),
      coef = c(1, -2, 1),
      u2 = c(0.1, 0.9),
      coef2 = c(-1, 1)
    )
  )

  expect_equal(res1$conf.int, res2$conf.int)
})


test_that("Groeneveld right skew versions agree", {

  res1 <- suppressWarnings(
    q.test(c(1:10), measure = "groenR")
  )

  res2 <- suppressWarnings(
    q.test(
      c(1:10),
      u = c(0.25, 0.5, 0.75),
      coef = c(1, -2, 1),
      u2 = c(0.25, 0.5),
      coef2 = c(-1, 1)
    )
  )

  expect_equal(res1$conf.int, res2$conf.int)
})


test_that("Groeneveld left skew versions agree", {

  res1 <- suppressWarnings(
    q.test(c(1:10), measure = "groenL")
  )

  res2 <- suppressWarnings(
    q.test(
      c(1:10),
      u = c(0.25, 0.5, 0.75),
      coef = c(1, -2, 1),
      u2 = c(0.5, 0.75),
      coef2 = c(1, -1)
    )
  )

  expect_equal(res1$conf.int, res2$conf.int)
})


test_that("Moors predefined and user-defined versions agree", {

  res1 <- suppressWarnings(
    q.test(c(1:10), measure = "moors")
  )

  res2 <- suppressWarnings(
    q.test(
      c(1:10),
      u = c(1/8, 3/8, 5/8, 7/8),
      coef = c(-1, 1, -1, 1),
      u2 = c(2/8, 6/8),
      coef2 = c(-1, 1)
    )
  )

  expect_equal(res1$conf.int, res2$conf.int)
})


test_that("left quantile weight versions agree", {

  res1 <- suppressWarnings(
    q.test(c(1:10), measure = "lqw")
  )

  res2 <- suppressWarnings(
    q.test(
      c(1:10),
      u = c(0.25/2, 0.25, 0.75/2),
      coef = c(1, -2, 1),
      u2 = c(0.25/2, 0.75/2),
      coef2 = c(-1, 1)
    )
  )

  expect_equal(res1$conf.int, res2$conf.int)
})


test_that("right quantile weight versions agree", {

  res1 <- suppressWarnings(
    q.test(c(1:10), measure = "rqw")
  )

  res2 <- suppressWarnings(
    q.test(
      c(1:10),
      u = c(1 - 0.75/2, 0.75, (1 + 0.75)/2),
      coef = c(1, -2, 1),
      u2 = c(1 - 0.75/2, (1 + 0.75)/2),
      coef2 = c(-1, 1)
    )
  )

  expect_equal(res1$conf.int, res2$conf.int)
})


test_that("qr1090 predefined and user-defined versions agree", {

  res1 <- suppressWarnings(
    q.test(c(1:10), measure = "qr1090")
  )

  res2 <- suppressWarnings(
    q.test(c(1:10), u = 0.1, u2 = 0.9)
  )

  expect_equal(res1$conf.int, res2$conf.int)
})


test_that("qr2575 predefined and user-defined versions agree", {

  res1 <- suppressWarnings(
    q.test(c(1:10), measure = "qr2575")
  )

  res2 <- suppressWarnings(
    q.test(c(1:10), u = 0.25, u2 = 0.75)
  )

  expect_equal(res1$conf.int, res2$conf.int)
})


test_that("q.test supports both variance methods", {

  set.seed(1234)
  x <- rlnorm(100)

  res.qor <- q.test(
    x,
    measure = "iqr",
    var.method = "qor"
  )

  res.density <- q.test(
    x,
    measure = "iqr",
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


test_that("q.test output is correct", {

  expect_equal(
    q.test(c(1:10))$method,
    "One sample test of the median"
  )

  expect_equal(
    q.test(c(1:10), c(1:10))$method,
    "Two sample test of the median"
  )

  expect_equal(
    q.test(c(1:10))$data.name,
    "c(1:10)"
  )

  expect_equal(
    as.numeric(q.test(c(1:10), c(1:10))$statistic),
    0
  )

  expect_equal(
    as.numeric(q.test(c(1:10), c(1:10))$p.value),
    1
  )

  expect_equal(
    q.test(c(1:10), alternative = "less")$alternative,
    "less"
  )

  expect_equal(
    as.numeric(q.test(c(1:10))$estimate),
    5.5
  )

  expect_equal(
    as.numeric(
      suppressWarnings(
        q.test(
          c(1:10),
          u = 0.25,
          coef = 1,
          u2 = 0.5,
          coef2 = 1,
          log.transf = TRUE
        )
      )$estimate
    ),
    as.numeric(
      log(
        quantile(c(1:10), 0.25, type = 8) /
          quantile(c(1:10), 0.5, type = 8)
      )
    )
  )

  expect_equal(
    as.numeric(
      suppressWarnings(
        q.test(
          c(1:10),
          u = 0.25,
          coef = 1,
          u2 = 0.5,
          coef2 = 1,
          log.transf = TRUE,
          back.transf = TRUE
        )
      )$estimate
    ),
    as.numeric(
      quantile(c(1:10), 0.25, type = 8) /
        quantile(c(1:10), 0.5, type = 8)
    )
  )

  expect_equal(
    as.numeric(q.test(c(1:10), true.q = 5)$null.value),
    5
  )
})
