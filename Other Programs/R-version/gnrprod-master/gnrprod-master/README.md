
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `gnrprod`

`gnrprod` implements the nonparametric identification of gross output
production functions outlined in Gandhi, Navarro, and Rivers (2020). The
main wrapper function `gnrprod` estimates the parameters of production
functions and productivity. The current version of `gnrprod` supports
only one flexible input.

## Installation

To install the latest development version of from Github, run the
following from R:

``` r
install.packages("devtools")
library(devtools)
devtools::install_github("davidjin0/gnrprod")
library(gnrprod)
```

## Usage

This package features three main functions: `gnrprod`, `gnrflex`, and
`gnriv`. `gnrprod` performs the entire production function estimation
routine. `gnrflex` performs the first stage of the estimation routine
which returns the flexible input elasticity. `gnriv` performs the second
stage of the estimation routine which returns the fixed input
elasticities and productivity.

An example of the use of `gnrprod`:

``` r
require(gnrprod)
#> Loading required package: gnrprod

# load Colombian plant-level data
data <- colombian

# estimate production function parameters and productivity
gnr_fit <- gnrprod(output = "RGO", fixed = c("L", "K"), flex = "RI",
                   share = "share", id = "id", time = "year", data = data)

# print results
gnr_fit
#> Gross Output Function:
#>   output:  "RGO" 
#>   fixed inputs:  c("L", "K") 
#>   flexible inputs:  "RI" 
#>   data: data
#> 
#> Estimates:
#>      L      K     RI 
#> 0.2332 0.1124 0.6793
summary(gnr_fit)
#> Gross Output Function:
#>   output: "RGO"
#>   fixed inputs:  c("L", "K") 
#>   flexible inputs: "RI"
#>   data: data
#> 
#> Total Productivity:
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#>   0.1416  12.5365  14.1472  14.8638  16.0179 539.9703 
#> 
#> Estimates:
#>    Estimate
#> L     0.233
#> K     0.112
#> RI    0.679
#> 
#> First-Stage Convergence: TRUE with 22 iterations. Sum of Squared Residuals: 315.1547 
#> Second-Stage Convergence: TRUE with 75 iterations. Value: 0.002222362
```

Alternatively, one can use `gnrflex` and `gnriv`:

``` r
# estimate flexible input elasticities (first stage)
gnr_fs <- gnrflex(output = "RGO", fixed = c("L", "K"), flex = "RI",
                  share = "share", id = "id", time = "year", data = data)

# print estimate
gnr_fs
#> Flexible Input Elasticity:
#>     RI 
#> 0.6793 
#> 
#> First-Stage Sum of Squared Residuals: 315.1547
#> Convergence: TRUE

# estimate fixed input elasticities and productivity (second stage)
gnr_ss <- gnriv(object = gnr_fs)

# print estimates
gnr_ss
#> Fixed Input Elasticity:
#>      L      K 
#> 0.2332 0.1124 
#> 
#> Total Productivity:
#>   productivity     
#>  Min.   :  0.1416  
#>  1st Qu.: 12.5365  
#>  Median : 14.1472  
#>  Mean   : 14.8638  
#>  3rd Qu.: 16.0179  
#>  Max.   :539.9703  
#> 
#> Second-Stage Objective Function Value: 0.002222362
#> Convergence:  TRUE
```

## References

Gandhi, Amit, Salvador Navarro, and David Rivers. 2020. “On the
Identification of Gross Output Production Functions.” *Journal of
Political Economy*, 128(8): 2973-3016. <https://doi.org/10.1086/707736>.
