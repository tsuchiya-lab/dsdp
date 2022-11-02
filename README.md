
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dsdp

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/dsdp)](https://CRAN.R-project.org/package=dsdp)
[![R-CMD-check](https://github.com/tsuchiya-lab/dsdp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tsuchiya-lab/dsdp/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
![GitHub](https://img.shields.io/github/license/tsuchiya-lab/dsdp)
<!-- badges: end -->

The goal of dsdp is to estimate probability density functions from a
data set using a maximum likelihood method. The models of density
functions in use are familiar Gaussian or exponential distributions with
polynomial correction terms. To find an optimal model, we adopt a grid
search for parameters of base functions and degrees of polynomials,
together with semidefinite programming for coefficients of polynomials,
and then model selection is done by Akaike Information Criterion.

## Installation

``` r
## Install from CRAN
install.packages(dsdp)
```

You can install the development version of dsdp like so:

``` r
## Install from github
devtools::install_github("tsuchiya-lab/dsdp")
```

To install from source codes, the user needs appropriate compiler
toolchain, for example, rtools in windows, to compile them, along with
devtools package.

## Usage

Please refer to the tutorials.
