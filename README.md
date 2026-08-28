
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rquest

<!-- badges: start -->

[![R-CMD-check](https://github.com/shenal-dkumara/rquest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/shenal-dkumara/rquest/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The `rquest` package provides convenient functionality for researchers
to carry out hypothesis tests and obtain confidence intervals for
measures based on quantiles. This includes single quantiles (e.g., the
median), linear combinations of quantiles (such as the interquartile
range), ratios of linear combinations commonly found in skewness and
kurtosis measures, and quantile-based inequality measures. Importantly,
the package also allows users to define their own measures based on
quantiles and carry out hypothesis tests and construct confidence
intervals for these measures.

The package also provides functionality for quantile density estimation,
computation of quantile optimality ratios (QOR), covariance estimation
for quantiles and ratios of linear combinations of quantiles, and
inference for robust versions of the coefficient of variation.

## Key functionality

`rquest` provides a flexible framework for quantile-based statistical
inference. Key functionality includes:

- **Quantile inference:** hypothesis tests and confidence intervals for
  individual quantiles, linear combinations of quantiles, and ratios of
  linear combinations of quantiles.

- **User-defined measures:** specification of custom quantile-based
  measures, allowing users to carry out hypothesis tests and construct
  confidence intervals beyond the measures built into the package.

- **One- and two-sample inference:** comparison of quantiles and
  quantile-based measures within a single population or between two
  independent populations.

- **Robust measures of variation:** inference for robust versions of the
  coefficient of variation based on the MAD and IQR.

- **Quantile-based inequality measures:** inference for built-in
  inequality measures as well as user-defined quantile-based inequality
  measures.

- **Quantile density estimation:** estimation of quantile densities
  using QOR-based or density-based methods.

- **Flexible QOR computation:** computation of quantile optimality
  ratios using the generalized lambda distribution, lognormal, normal,
  exponential, or a user-supplied quantile function.

- **Covariance estimation:** covariance matrices for quantile estimators
  and for ratios of linear combinations of quantiles.

## Main functions

The main functions in the package are:

- `q.test()` carries out hypothesis tests and obtains associated
  confidence intervals for quantiles, linear combinations of quantiles,
  and ratios of linear combinations of quantiles.

- `qineq()` carries out hypothesis tests and obtains associated
  confidence intervals for quantile-based inequality measures, including
  user-defined measures.

- `qcov()` computes covariance matrices containing variances and
  covariances of quantile estimators.

- `qden()` estimates the quantile density function using either a
  QOR-based method or a density-based method.

- `qor()` computes quantile optimality ratios for supported
  distributions or user-supplied quantile functions.

- `qrcov()` computes approximate covariance matrices for ratios of
  linear combinations of quantile estimators.

- `rcv.test()` carries out one- and two-sample inference for robust
  versions of the coefficient of variation based on the MAD or IQR.

## Installation

You can install the released version of `rquest` from CRAN with:

``` r
install.packages("rquest")
```

You can install the development version of `rquest` from GitHub with:

``` r
# install.packages("pak")
pak::pak("shenal-dkumara/rquest")
```

## Usage

``` r
library(rquest)
```

### Quantile-based inference with `q.test()`

``` r
# Create some data
x <- c(8.43, 7.08, 8.79, 8.88, 7.87,
       5.94, 8.79, 5.46, 8.11, 7.08)

y <- c(13.44, 13.65, 14.77, 9.51, 14.07,
       10.92, 11.59, 13.42, 8.93, 10.88)

# One-sample hypothesis test for the IQR
q.test(x, measure = "iqr")
#> 
#>  One sample test of the interquartile range (IQR)
#> 
#> data:  x
#> Z = 1.4037, p-value = 0.1604
#> alternative hypothesis: true IQR  is not equal to 0
#> 95 percent confidence interval:
#>  -0.7153314  4.3253314
#> sample estimates:
#>  IQR  
#> 1.805

# Two-sample hypothesis test for robust coefficients of variation
# (0.75 * IQR / median), with log transformation and
# back-transformation to the ratio scale
q.test(
  x, y,
  measure = "rCViqr",
  log.transf = TRUE,
  back.transf = TRUE
)
#> 
#>  Two sample test of the robust coefficient of variation
#>  (0.75*IQR/median)
#> 
#> data:  x and y
#> Z = -0.032603, p-value = 0.974
#> alternative hypothesis: true ratio of Robust CVs  is not equal to 1
#> 95 percent confidence interval:
#>  0.1347132 6.9518490
#> sample estimates:
#> ratio of Robust CVs  
#>            0.9677323
```

The same measure can also be specified directly in terms of the required
quantiles and coefficients:

``` r
u <- c(0.25, 0.75)
coef <- 0.75 * c(-1, 1)
u2 <- 0.5
coef2 <- 1

q.test(
  x, y,
  u = u,
  u2 = u2,
  coef = coef,
  coef2 = coef2,
  log.transf = TRUE,
  back.transf = TRUE
)
#> 
#>  Two sample test of a ratio of two linear combinations of quantiles
#>  (LCQs)
#> 
#> data:  x and y
#> Z = -0.032603, p-value = 0.974
#> alternative hypothesis: true ratio of Ratio of LCQs  is not equal to 1
#> 95 percent confidence interval:
#>  0.1347132 6.9518490
#> sample estimates:
#> ratio of Ratio of LCQs  
#>               0.9677323
```

Alternatively, the numerator and denominator can be defined using a
coefficient matrix:

``` r
u <- c(0.25, 0.5, 0.75)
num <- 0.75 * c(-1, 0, 1)
den <- c(0, 1, 0)
coef <- rbind(num, den)

q.test(
  x, y,
  u = u,
  coef = coef,
  log.transf = TRUE,
  back.transf = TRUE
)
#> 
#>  Two sample test of a ratio of two linear combinations of quantiles
#>  (LCQs)
#> 
#> data:  x and y
#> Z = -0.032603, p-value = 0.974
#> alternative hypothesis: true ratio of Ratio of LCQs  is not equal to 1
#> 95 percent confidence interval:
#>  0.1347132 6.9518490
#> sample estimates:
#> ratio of Ratio of LCQs  
#>               0.9677323
```

### Covariance estimation with `qcov()`

``` r
set.seed(1234)
xq <- rnorm(100)

# Covariance matrix for sample quartiles
qcov(xq, c(0.25, 0.5, 0.75))
#>             0.25         0.5        0.75
#> 0.25 0.014141467 0.008979251 0.008130512
#> 0.5  0.008979251 0.017104367 0.015487625
#> 0.75 0.008130512 0.015487625 0.042071102
```

### Quantile density estimation with `qden()`

``` r
set.seed(1234)
xd <- rnorm(100)

# QOR-based quantile density estimation using
# the flexible GLD distribution
qden(
  xd,
  c(0.25, 0.5, 0.75),
  method = "qor"
)
#> [1] 2.746291 2.615673 4.736868

# QOR-based quantile density estimation using
# the normal distribution
qden(
  xd,
  c(0.25, 0.5, 0.75),
  dist = "norm",
  method = "qor"
)
#> [1] 2.409224 2.638314 4.776829

# Density-based quantile density estimation
qden(
  xd,
  c(0.25, 0.5, 0.75),
  method = "density"
)
#> [1] 2.542878 2.386256 4.299268
```

### Quantile optimality ratios with `qor()`

``` r
set.seed(1234)
xqor <- rlnorm(100)

# QOR using the flexible FKML GLD distribution
# with parameters estimated from the data
qor(
  seq(0.1, 0.9, by = 0.2),
  x = xqor
)
#> $qor
#> [1] -0.458044160  0.214678367  0.088987175  0.027832449  0.002757499
#> 
#> $params
#> NULL

# QOR using the lognormal distribution
qor(
  seq(0.1, 0.9, by = 0.2),
  dist = "lnorm",
  x = xqor
)
#> $qor
#> [1] 0.021508337 0.123522437 0.079227681 0.029207901 0.003364183
#> 
#> $params
#> NULL

# QOR for the exponential distribution
qor(
  seq(0.1, 0.9, by = 0.2),
  dist = "exp"
)
#> $qor
#> [1] 0.405 0.245 0.125 0.045 0.005
#> 
#> $params
#> NULL

# QOR for a user-supplied quantile function
qor(
  c(0.2, 0.5, 0.8),
  dist = qbeta,
  params = list(shape1 = 2, shape2 = 5)
)
#> $qor
#> [1] 0.07645133 0.19841722 0.03361339
#> 
#> $params
#> $params$shape1
#> [1] 2
#> 
#> $params$shape2
#> [1] 5
```

### Covariance estimation for ratios with `qrcov()`

``` r
set.seed(1234)
xr <- rnorm(100)

# Covariance matrix for the ratio of the third quartile
# to the median and the ratio of the first quartile
# to the median
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

qrcov(
  xr,
  c(0.25, 0.5, 0.75),
  coef1 = coef1,
  coef2 = coef2
)
#> $ratios
#> [1] -1.265365  2.346111
#> 
#> $cov
#>            R1         R2
#> R1  0.7344446 -0.4570857
#> R2 -0.4570857  0.4471806
```

### Quantile-based inequality measures with `qineq()`

``` r
# Create some data
x <- c(8.43, 7.08, 8.79, 8.88, 7.87,
       5.94, 8.79, 5.46, 8.11, 7.08)

y <- c(13.44, 13.65, 14.77, 9.51, 14.07,
       10.92, 11.59, 13.42, 8.93, 10.88)

# Test equality of the QRI measure between two groups
qineq(x, y, measure = "QRI")
#> 
#>  Two sample test of the QRI statistic
#> 
#> data:  x and y
#> Z = -0.21439, p-value = 0.8302
#> alternative hypothesis: true difference in QRI statistic is not equal to 0
#> 95 percent confidence interval:
#>  -0.2522859  0.2025363
#> sample estimates:
#> difference in QRI statistic 
#>                 -0.02487481
```

A user-defined inequality measure can be specified directly using a
formula. For example, the QRI can equivalently be defined as:

``` r
qineq(
  x, y,
  measure = ~ Q(p/2) / Q(1 - p/2)
)
#> 
#>  Two sample test of the user defined statistic
#> 
#> data:  x and y
#> Z = -0.21439, p-value = 0.8302
#> alternative hypothesis: true difference in user defined statistic is not equal to 0
#> 95 percent confidence interval:
#>  -0.2522859  0.2025363
#> sample estimates:
#> difference in user defined statistic 
#>                          -0.02487481
```

More generally, a user-defined measure can also be specified directly
using the required quantiles and coefficient matrices:

``` r
J <- 100
p <- (1:J - 0.5) / J
u <- sort(c(p/2, 1 - p/2))

num <- cbind(
  diag(rep(1, J)),
  matrix(0, J, J)
)

den <- cbind(
  matrix(0, J, J),
  diag(J)[, J:1]
)

qineq(
  x, y,
  measure = list(
    u = u,
    coef1 = num,
    coef2 = den
  )
)
#> 
#>  Two sample test of the user defined statistic
#> 
#> data:  x and y
#> Z = -0.21439, p-value = 0.8302
#> alternative hypothesis: true difference in user defined statistic is not equal to 0
#> 95 percent confidence interval:
#>  -0.2522859  0.2025363
#> sample estimates:
#> difference in user defined statistic 
#>                          -0.02487481
```

### Robust coefficients of variation with `rcv.test()`

``` r
set.seed(123)
x1 <- rnorm(100, 8)
x2 <- rnorm(120, 10)

# One-sample inference using the MAD
rcv.test(x1)
#> 
#>  One sample test of the robust coefficient of variation (MAD/median)
#> 
#> data:  x1
#> Z = -18.082, p-value < 2.2e-16
#> alternative hypothesis: true Robust CV  is not equal to 1
#> 95 percent confidence interval:
#>  0.08691013 0.14014523
#> sample estimates:
#> Robust CV  
#>  0.1103632

# One-sample inference using the IQR
rcv.test(x1, numerator = "iqr")
#> 
#>  One sample test of the robust coefficient of variation
#>  (0.75*IQR/median)
#> 
#> data:  x
#> Z = -16.613, p-value < 2.2e-16
#> alternative hypothesis: true Robust CV  is not equal to 1
#> 95 percent confidence interval:
#>  0.0856860 0.1439217
#> sample estimates:
#> Robust CV  
#>  0.1110499

# Two-sample inference using the MAD
rcv.test(x1, x2)
#> 
#>  Two sample test of the robust coefficient of variation (MAD/median)
#> 
#> data:  x1 and x2
#> Z = 0.97732, p-value = 0.3284
#> alternative hypothesis: true ratio of Robust CVs  is not equal to 1
#> 95 percent confidence interval:
#>  0.8489523 1.6314679
#> sample estimates:
#> ratio of Robust CVs  
#>             1.176877

# Two-sample inference using the IQR
rcv.test(x1, x2, numerator = "iqr")
#> 
#>  Two sample test of the robust coefficient of variation
#>  (0.75*IQR/median)
#> 
#> data:  x and y
#> Z = 0.65049, p-value = 0.5154
#> alternative hypothesis: true ratio of Robust CVs  is not equal to 1
#> 95 percent confidence interval:
#>  0.7898182 1.6005955
#> sample estimates:
#> ratio of Robust CVs  
#>             1.124357
```
