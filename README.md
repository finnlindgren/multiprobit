
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multiprobit

<!-- badges: start -->

[![Travis build status of devel
branch](https://travis-ci.org/finnlindgren/multiprobit.svg?branch=devel)](https://travis-ci.org/finnlindgren/multiprobit)
[![Codecov test
coverage](https://codecov.io/gh/finnlindgren/multiprobit/branch/master/graph/badge.svg)](https://codecov.io/gh/finnlindgren/multiprobit?branch=master)
<!-- badges: end -->

The goal of multiprobit is to perform fast Bayesian inference for
multivariate probit models. The method uses a latent Gaussian variable
parameterisation of the correlation matrix, and numerical optimisation
and integration to find the posterior distributions of the model
coefficients and correlation matrix.

## Installation

<!--
You can install the released version of `multiprobit` from [CRAN](https://CRAN.R-project.org) with:


```r
install.packages("multiprobit")
```
-->

To install the development version of `multiprobit` from github, use

``` r
remotes::install_github("finnlindgren/multiprobit", ref = "devel")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(multiprobit)
## basic example code
```
