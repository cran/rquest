
test_that("type", {
  expect_true(is.matrix(qcov(c(1:10),c(0.25, 0.5, 0.75))))
})


test_that("variance", {
  expect_true(
  all.equal(qcov(c(1:10), c(0.25, 0.25))[1,1],
            qcov(c(1:10), c(0.25, 0.25))[1,2],
            qcov(c(1:10), c(0.25, 0.25))[2,1],
            qcov(c(1:10), c(0.25, 0.25))[2,2],
            qcov(c(1:10), c(0.25))[1,1]))
  expect_equal(as.numeric(qcov(rep(1,10),0.25)), 0)
})







