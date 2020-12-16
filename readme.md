
# echoice2 <a href='http://ninohardt.de/echoice2/'><img src="c:/users/myself/RGithub/echoice2/man/figures/echoicelogo.png" align="right" height="139"/></a>

## Overview

This package contains choice models with economic foundation.

For more theoretical background, please refer to the chapter ‘[Economic
foundations of conjoint
analysis](https://doi.org/10.1016/bs.hem.2019.04.002)’ in the [Handbook
of the Economics of
Marketing](https://www.elsevier.com/books/handbook-of-the-economics-of-marketing/dube/978-0-444-63759-8).

For details on how to use echoice, please refer to the
[vignette](http://ninohardt.de/echoice2/articles/echoice2.html) in the
package. It illustrates model estimation, predictions and model
evaluation.

All key functions are written in c++ and use
[openMP](https://www.openmp.org) for multi-threaded computing. C++
integration in R is facilitated by [Rcpp](http://www.rcpp.org).

`echoice2` follows tidy principles and integrated nicely with `dplyr`.
It can be used to generate choice volume/share simulators, though no
front-end is built into the package yet.

`echoice2` is a complete rewrite of the original `echoice` package,
which is still available on
[Bitbucket](https://bitbucket.org/ninohardt/echoice).

## Installation

``` r
# Install the cutting edge development version from GitHub:
# install.packages("devtools")
devtools::install_github("ninohardt/echoice2")
#remotes::install_github("ninohardt/echoice2")
```

### Installation notes

-   Make sure you are using a recent R4.0.x version, `tidyerse`, `Rcpp`
    and `RcppArmadillo` before building the package from Github.

-   If you use Linux, this should compile just fine.

-   If you are using OSX, you may have to install CLI, XQuartz and
    potentially other things that Apple removed from OSX. Google ‘+OSX
    +Rcpp’ if you run into trouble.

-   If you are using Windows, install
    \[Rtools\](<https://cran.r-project.org/bin/windows/Rtools/)> first.

### Binaries

Binaries for windows can be downloaded
[here](http://ninohardt.de/echoice2/echoice2_0.1.zip "echoice2"). You
can install the binary version from the ‘Packages’ tab in RStudio.
Select ‘Install’ and change ‘Install from’ to ‘Package Archive File’. Or
you use \`install.packages’ on the command line and point it to the
downloaded .zip file.

## Functionality

The following models are implemented (including estimation and
prediction):

-   Discrete Choice (HMNL)
-   Discrete Choice, attribute-based screening (not including price)
-   Volumetric demand, EV1 errors
-   Volumetric demand, Normal errors
-   Volumetric demand, attribute-based screening, Normal errors
-   Volumetric demand, attribute-based screening including price, Normal
    errors
-   Volumetric demand, accounting for set-size variation (1st order),
    EV1 errors
-   Volumetric demand, accounting for set-size variation (2nd order),
    EV1 errors

## Usage

Please read the vignette. It illustrates a complete workflow from
estimation to market simulation. The vignette can also be found on the
[package website](http://ninohardt.de/echoice2/articles/echoice2.html).

Functions that relate to discrete demand start in `dd_`, while functions
for volumetric demand start in `vd_`. Universal functions (discrete and
volumetric choice) start in `ec_`. Estimation functions continue in
`est`, demand simulators in `dem.`

The package comes with a small example dataset `icecream` from a
volumetric conjoint study. It contains 300 respondents.

``` r
  data(icecream)
  icecream %>% head
#> # A tibble: 6 x 8
#>      id  task   alt     x     p Brand     Flavor      Size 
#>   <int> <int> <int> <dbl> <dbl> <fct>     <fct>       <ord>
#> 1     1     1     1     8 0.998 Store     Neapolitan  16   
#> 2     1     1     2     0 0.748 Store     VanillaBean 16   
#> 3     1     1     3     0 1.25  BenNJerry Oreo        16   
#> 4     1     1     4     0 0.748 BenNJerry Neapolitan  16   
#> 5     1     1     5     0 2.49  HaagenDa  RockyRoad   4    
#> 6     1     1     6     0 1.25  HaagenDa  Oreo        16
```

Choice data data.frames or tibbles need to contain the following
columns:

-   `id` (integer; respondent identifier)
-   `task` (integer; task number)
-   `alt` (integer; alternative number within task)
-   `x` (double; quantity purchased)
-   `p` (double; price)
-   attributes defining the choice alternatives (factor)

While this requires a little extra space for discrete choice data, it
simplifies the workflow.

Estimating a simple volumetric demand model is easy. Use the
`vd_est_vdm` function:

``` r
icecream <- ice_cal %>% vd_est_vdm(R=10000)
```
