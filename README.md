
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multiprobit

<!-- badges: start -->

[![R build
status](https://github.com/finnlindgren/multiprobit/workflows/R-CMD-check/badge.svg)](https://github.com/finnlindgren/multiprobit/actions)
[![R code coverage
status](https://github.com/finnlindgren/multiprobit/workflows/test-coverage/badge.svg)](https://github.com/finnlindgren/multiprobit/actions)
[![Codecov test
coverage](https://codecov.io/gh/finnlindgren/multiprobit/branch/devel/graph/badge.svg)](https://codecov.io/gh/finnlindgren/multiprobit?branch=devel)
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
if (interactive()) {
  library(multiprobit)
  
  N <- 6
  d <- 2
  J <- 2
  
  set.seed(1L)
  X <- cbind(1, matrix(rnorm(N * (J - 1)), N, J - 1))
  B <- matrix(0.5, J, d)
  Y <- matrix(rnorm(N * d, mean = as.vector(X %*% B)) > 0, N, d)
  df <- d + 1
  prec_beta <- 0.1
  
  model <- mp_model(
    response = Y, X = X,
    df = df, prec_beta = prec_beta
  )
  opt <- multiprobit(
    model = model,
    options =
      mp_options(
        gaussint = list(max.threads = 1),
        strategy = "stepwise"
      )
  )
}
```
