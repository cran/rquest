test_that("class", {
  expect_s3_class(qineq(c(1:10)), "htest")
})

test_that("remove NA values: warning", {
  expect_warning(qineq(c(1, 2, NA, 4)))
  expect_warning(qineq(c(1, 2, NA, 4),c(1, 2, NA, 4)))
})

test_that("errors in arguments", {
  expect_error(qineq(c(1:10),measure = "Gini"),"Only QRI and the G2 measures are supported at the moment.")
})

test_that("output", {
  expect_equal(qineq(c(1:10))$method,"One sample test of the QRI")
  expect_equal(qineq(c(1:10),c(1:10))$method,"Two sample test of the QRI")
  expect_equal(qineq(c(1:10),measure = "G2")$method,"One sample test of the G2")
  expect_equal(qineq(c(1:10))$data.name,"c(1:10)")
  expect_equal(as.numeric(qineq(c(1:10),c(1:10))$statistic),0)
  expect_equal(qineq(c(1:10),alternative = "less")$alternative,"less")
  expect_equal(as.numeric(qineq(c(1:10),c(1:10))$estimate),0)
  expect_equal(as.numeric(qineq(c(rep(1,10)))$estimate),0)
  expect_equal(as.numeric(qineq(c(1:10),c(1:10),true.ineq = 5)$null.value),5)
})
