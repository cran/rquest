test_that("class", {
  expect_s3_class(q.test(c(1:10)), "htest")
})

test_that("remove NA values: warning", {
  expect_warning( q.test(c(1, 2, NA, 4)))
  expect_warning(q.test(c(1, 2, NA, 4),c(1, 2, NA, 4)))
  expect_warning(q.test(c(1:10),u=c(0.25,0.5),coef=c(-1,1),u2=c(0.5,0.75),coef2=c(-1,1)))
  expect_no_warning(q.test(c(1:10),u=c(0.25,0.5),coef=c(-1,1),u2=c(0.5,0.75),coef2=c(-1,1),log.transf = TRUE))
})

test_that("log transformation: warning", {
  expect_warning(q.test(c(1, 2, NA, 4)))
  expect_warning(q.test(c(1, 2, NA, 4),c(1, 2, NA, 4)))
})


test_that("errors in arguments", {
  expect_error(q.test(c(1:10),coef2=1),"When using u2 and coef2, you also need to specify both u and coef.")
  expect_error(q.test(c(1:10),coef=1,coef2=1),"When using u2 and coef2, you also need to specify both u and coef.")
  expect_error(q.test(c(1:10),measure = "bowley",p=0.5),"Argument p must be a numeric value in (0, 1) except 0.5.",fixed=TRUE)
  expect_error(q.test(c(1:10),measure="measure"),"Unknown choice for measure.")
  expect_error(q.test(c(1:10),measure = "groen"),"Use one of 'groenR' or 'groenL'.")
  expect_error(q.test(c(1:10),measure = "groenL",p=0.5),"Argument p must be a numeric value in (0, 1) except 0.5.",fixed=TRUE)
  expect_error(q.test(c(1:10),measure = "lqw",p=0.75),"Argument p must be a numeric value in (0, 1/2).",fixed=TRUE)
  expect_error(q.test(c(1:10),measure = "rqw",p=0.25),"Argument p must be a numeric value in (1/2, 1).",fixed=TRUE)
  expect_error(q.test(c(1:10),measure="qr0.250.75"),"For quantile ratios, measure must be in format qrxxyy where xx and yy are integer numbers.")
  expect_error(q.test(c(1:10),u=c(0.25,0.75),coef = c(1,"a")),"If coef is a vector, then it needs to be a numeric vector.")
  expect_error(q.test(c(1:10),coef=1),"Argument u required if coef is a numeric vector.")
  expect_error(q.test(c(1:10),u=0.25,coef=c(1,-1)),"Length of u needs to be equal to the length of coef.")
  expect_error(q.test(c(1:10),u=c(0.25,0.75),coef = matrix(c(1,"a",0,1),2,2)),"If coef is a matrix, then it needs to be a numeric matrix.")
  expect_error(q.test(c(1:10),coef = matrix(c(1,0,0,1),2,2)),"Argument u required if coef is a numeric matrix.")
  expect_error(q.test(c(1:10),u=c(0.25,0.75),coef = matrix(c(1,0,1,0,1,1),3,2)),"Matrix coef needs to have dimensions: ncol(coef)=length(u) and nrow(coef)=2.",fixed=TRUE)
  expect_error(q.test(c(1:10),var.method = "method"),"Argument method must be either 'qor' or 'density'.")
})

test_that("check LCQ is correct", {
  expect_equal(as.numeric(q.test(c(1:10),u=c(0.25,0.75),coef=c(-2,4),u2=c(0.5,0.75),coef2=c(-2,3))$estimate),2)
})

test_that("vector and matrix inputs", {
  expect_equal(q.test(c(1:10),u=c(0.25,0.75),coef=c(-1,1),u2=c(0.5,0.75),coef2=c(-1,1)),
               q.test(c(1:10),u=c(0.25,0.5,0.75),coef=matrix(c(-1,0,1,0,-1,1), 2, 3,byrow=TRUE)))
})

test_that("rCViqr", {
  expect_equal(q.test(c(1:10),measure = "rCViqr")$conf.int,
               q.test(c(1:10),u=c(0.25,0.75),coef=c(-0.75,0.75),u2=0.5)$conf.int)
})

test_that("IQR", {
  expect_equal(q.test(c(1:10),measure = "IQR")$conf.int,
               q.test(c(1:10),u=c(0.25,0.75),coef=c(-1,1))$conf.int)
})

test_that("Median", {
  expect_equal(q.test(c(1:10),measure = "median")$conf.int,
               q.test(c(1:10),u=0.5)$conf.int)
})

test_that("bowley", {
  expect_equal(q.test(c(1:10),measure = "bowley")$conf.int,
               q.test(c(1:10),u=c(0.25,0.5,0.75),coef=c(1,-2,1),u2=c(0.25,0.75),coef2=c(-1,1))$conf.int)
})

test_that("kelly", {
  expect_equal(q.test(c(1:10),measure = "kelly")$conf.int,
               q.test(c(1:10),u=c(0.1,0.5,0.9),coef=c(1,-2,1),u2=c(0.1,0.9),coef2=c(-1,1))$conf.int)
})


test_that("groenR", {
  expect_equal(q.test(c(1:10),measure = "groenR")$conf.int,
               q.test(c(1:10),u=c(0.25,0.5,0.75),coef=c(1,-2,1),u2=c(0.25,0.5),coef2=c(-1,1))$conf.int)
})

test_that("groenL", {
  expect_equal(q.test(c(1:10),measure = "groenL")$conf.int,
               q.test(c(1:10),u=c(0.25,0.5,0.75),coef=c(1,-2,1),u2=c(0.5,0.75),coef2=c(1,-1))$conf.int)
})



test_that("moors", {
  expect_equal(q.test(c(1:10),measure = "lqw")$conf.int,
               q.test(c(1:10),u=c(0.25/2,0.25,0.75/2),coef=c(1,-2,1),u2=c(0.25/2,0.75/2),coef2=c(-1,1))$conf.int)
})


test_that("lqw", {
  expect_equal(q.test(c(1:10),measure = "moors")$conf.int,
               q.test(c(1:10),u=c(1/8,3/8,5/8,7/8),coef=c(-1,1,-1,1),u2=c(2/8,6/8),coef2=c(-1,1))$conf.int)
})

test_that("rqw", {
  expect_equal(q.test(c(1:10),measure = "rqw")$conf.int,
               q.test(c(1:10),u=c(1-0.75/2,0.75,(1+0.75)/2),coef=c(1,-2,1),u2=c(1-0.75/2,(1+0.75)/2),coef2=c(-1,1))$conf.int)
})


test_that("qr1090", {
  expect_equal(q.test(c(1:10),measure = "qr1090")$conf.int,
               q.test(c(1:10),u=0.1,u2=0.9)$conf.int)
})

test_that("equivalence", {
  expect_equal(q.test(c(1:10),measure = "qr2575")$conf.int,
               q.test(c(1:10),u=0.25,u2=0.75)$conf.int)
})


test_that("output", {
  expect_equal(q.test(c(1:10))$method,"One sample test of the median")
  expect_equal(q.test(c(1:10),c(1:10))$method,"Two sample test of the median")
  expect_equal(q.test(c(1:10))$data.name,"c(1:10)")
  expect_equal(as.numeric(q.test(c(1:10),c(1:10))$statistic),0)
  expect_equal(as.numeric(q.test(c(1:10),c(1:10))$p.value),1)
  expect_equal(q.test(c(1:10),alternative = "less")$alternative,"less")
  expect_equal(as.numeric(q.test(c(1:10))$estimate),5.5)
  expect_equal(as.numeric(q.test(c(1:10),u=0.25,coef=1,u2=0.5,coef2=1,log.transf = TRUE)$estimate),as.numeric(log(quantile(c(1:10),0.25,type = 8)/quantile(c(1:10),0.5,type = 8))))
  expect_equal(as.numeric(q.test(c(1:10),u=0.25,coef=1,u2=0.5,coef2=1,log.transf = TRUE,back.transf = TRUE)$estimate),as.numeric(quantile(c(1:10),0.25,type = 8)/quantile(c(1:10),0.5,type = 8)))
  expect_equal(as.numeric(q.test(c(1:10),true.q = 5)$null.value),5)

})





