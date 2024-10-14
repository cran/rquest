
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rquest

<!-- badges: start -->

[![R-CMD-check](https://github.com/shenal-dkumara/rquest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/shenal-dkumara/rquest/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The `rquest` package provides convenient functionality for researchers
to carry out hypothesis tests and obtain confidence intervals for
measures based on quantiles. This includes for single quantiles (e.g.,
the median), linear combinations of quantiles (such as the interquartile
range), ratios of linear combinations commonly found in skewness and
kurtosis measures, and newly developed inequality measures. Another key
objective is to make it easy for users to define their own measures for
hypothesis testing and confidence intervals.

Following are the main functions in the package:

- `q.test()` carry out hypothesis tests and obtain associated confidence
  intervals for linear combinations of quantiles, and ratios of such
  linear combinations.

- `qineq()` carry out hypothesis tests and obtain associated confidence
  intervals for quantile based inequality measures

- `qcov()` compute a covariance matrix consisting of variances (on the
  diagonal) for quantile estimates and covariances (off-diagonal)
  between different quantile estimates

## Installation

You can install the development version of `rquest` from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("shenal-dkumara/rquest")
```

## Usage

``` r
library(rquest)

## Functionality of q.test() ##

#  Create some data
x <- c(8.43,7.08,8.79,8.88,7.87,5.94,8.79,5.46,8.11,7.08)
y <- c(13.44,13.65,14.77,9.51,14.07,10.92,11.59,13.42,8.93,10.88)

# One sample hypothesis test for the IQR
q.test(x, measure = "iqr")
#> 
#>  One sample test of the interquartile range (IQR)
#> 
#> data:  x
#> Z = 2.4436, p-value = 0.01454
#> alternative hypothesis: true IQR  is not equal to 0
#> 95 percent confidence interval:
#>  0.3572545 3.2527455
#> sample estimates:
#>  IQR  
#> 1.805

# Two samples hypothesis test for robust coefficient variations (0.75*IQR/median) with log transformation and back-transformation to the ratio scale,.
q.test(x, y, measure = "rCViqr", log.transf = TRUE, back.transf = TRUE)
#> 
#>  Two sample test of the robust coefficient of variation
#>  (0.75*IQR/median)
#> 
#> data:  x and y
#> Z = -0.059465, p-value = 0.9526
#> alternative hypothesis: true ratio of Robust CVs  is not equal to 1
#> 95 percent confidence interval:
#>  0.3282838 2.8527321
#> sample estimates:
#> ratio of Robust CVs  
#>            0.9677323

# The same two samples hypothesis test for robust coefficient variations (0.75*IQR/median) by using 'u',''u2','coef' and 'coef2' arguments.
u<-c(0.25,0.75)
coef<-0.75*c(-1,1)
u2<-0.5
coef2<-1
q.test(x,y,u=u,u2=u2,coef=coef,coef2=coef2,log.transf=TRUE,back.transf=TRUE)
#> 
#>  Two sample test of a ratio of two linear combinations of quantiles
#>  (LCQs)
#> 
#> data:  x and y
#> Z = -0.059465, p-value = 0.9526
#> alternative hypothesis: true ratio of Ratio of LCQs  is not equal to 1
#> 95 percent confidence interval:
#>  0.3282838 2.8527321
#> sample estimates:
#> ratio of Ratio of LCQs  
#>               0.9677323

# The same two samples hypothesis test for robust coefficient variations (0.75*IQR/median) by using only 'u' and 'coef' arguments.
u<-c(0.25,0.5,0.75)
num <- 0.75*c(-1,0,1)
den <- c(0,1,0)
coef <- rbind(num, den)
q.test(x,y,u=u,coef=coef,log.transf=TRUE,back.transf=TRUE)
#> 
#>  Two sample test of a ratio of two linear combinations of quantiles
#>  (LCQs)
#> 
#> data:  x and y
#> Z = -0.059465, p-value = 0.9526
#> alternative hypothesis: true ratio of Ratio of LCQs  is not equal to 1
#> 95 percent confidence interval:
#>  0.3282838 2.8527321
#> sample estimates:
#> ratio of Ratio of LCQs  
#>               0.9677323

## Functionality of qcov() ##

# Compute the variance-covariance matrix for sample quartiles.
qcov(x, c(0.25, 0.5, 0.75))
#> $cov
#>            0.25        0.5       0.75
#> 0.25 0.61381325 0.27565677 0.06757741
#> 0.5  0.27565677 0.37138326 0.09104481
#> 0.75 0.06757741 0.09104481 0.06695905


## Functionality of qineq() ##

# Two sample hypothesis test for the QRI measure
qineq(x,y)
#> 
#>  Two sample test of the QRI
#> 
#> data:  x and y
#> Z = -0.28062, p-value = 0.779
#> alternative hypothesis: true difference in QRI is not equal to 0
#> 95 percent confidence interval:
#>  -0.2158219  0.1617608
#> sample estimates:
#> difference in QRI 
#>       -0.02703056
```
