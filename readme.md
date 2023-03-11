
# echoice2 <a href='http://ninohardt.de/echoice2/'><img src="http://ninohardt.de/echoice2/reference/figures/echoicelogo.png" align="right" height="139"/></a>

## Overview

This package contains choice models with economic foundation. Its
purpose is to simplify using choice models with economic foundation. Key
tenets are: (1) Simple, flexible data handling that is compatible with
R-tidyverse and general enough to support many different models (2)
speed.

For more theoretical background and reasons to use choice models with
economic foundation, please refer to the chapter ‘[Economic foundations
of conjoint analysis](https://doi.org/10.1016/bs.hem.2019.04.002)’ in
the [Handbook of the Economics of
Marketing](https://www.elsevier.com/books/handbook-of-the-economics-of-marketing/dube/978-0-444-63759-8).

All key functions are written in c++ and use
[openMP](https://www.openmp.org) for multi-threaded computing. C++
integration in R is facilitated by [Rcpp](http://www.rcpp.org).

`echoice2` (largely) follows tidy principles and integrated nicely with
`dplyr`. It can be used to generate choice volume/share simulators,
though no front-end is built into the package yet.

## News

- **Looking for co-developers!**

- Version
  [**0.22**](https://github.com/ninohardt/echoice2/releases/tag/v0.2.2):

  - Added a vignette on volumetric choice modeling with and without 
    conjunctive screening
  - CRAN release to follow soon
    
- Version
  [**0.21**](https://github.com/ninohardt/echoice2/releases/tag/v0.2.1):

  - Package compiles in absense of OpenMP Support
  - MacOS CRAN binaries are compiled without OpenMP - for full speed
    compile yourself or use binaries supplied here

- Version
  [**0.20**](https://github.com/ninohardt/echoice2/releases/tag/v0.20):

  - Initial CRAN release
  - No new functionality, but cleaner code

- Version
  [**0.16**](https://github.com/ninohardt/echoice2/releases/tag/v0.16):

  - Stability/Performance improvements for demand simulators

- Version
  [**0.15**](https://github.com/ninohardt/echoice2/releases/tag/v0.15):

  - Initial release
  - faster and more efficient screening model estimation
  - improved demand predictions: posterior demand draws are now stored
    in a single column - this is a *major* improvement for dealing with
    demand predictions!
  - some bug-fixes and documentation improvements

## Installation

``` r
#install from CRAN
# install.packages("echoice2")

#install from github
# install.packages("remotes")
# remotes::install_github("ninohardt/echoice2", build_vignettes = TRUE)
```

``` r
library(echoice2)
```

### Installation notes

#### If installing from source/github

- If you use Linux it should just work.
- If you are using Windows, install
  [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first.
- If you are using OSX, you may have to install CLI, XQuartz and
  potentially other things that Apple removed from OSX. For multicore
  support, you also need to find a compiler that does support OpenMP.
  Just google it.

## Functionality

The following models are implemented (including estimation and
prediction):

- Discrete Choice (HMNL)
  - Without Screening
  - With conjunctive Screening
- Volumetric demand (EV1, Normal errors)
  - Without Screening
  - With conjunctive Screening
  - With set-size variation

## Ideas about future functionality

- upper-level covariates, effects-coding, discrete and continuous
  attributes, discrete choice example vignette, …

## Usage

Functions that relate to discrete demand start in `dd_`, while functions
for volumetric demand start in `vd_`. Universal functions (discrete and
volumetric choice) start in `ec_`. Estimation functions continue in
`est`, demand simulators in `dem.`

The package comes with a small example dataset `icecream` from a
volumetric conjoint study. It contains 300 respondents.

``` r
  data(icecream)
  icecream %>% head
#> # A tibble: 6 × 8
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

- `id` (integer; respondent identifier)
- `task` (integer; task number)
- `alt` (integer; alternative number within task)
- `x` (double; quantity purchased)
- `p` (double; price)
- attributes defining the choice alternatives (factor, and soon
  continuous as well)

While this requires a little extra space for discrete choice data, it
simplifies the workflow and makes the package versatile. It can be
applied to data from choice experiments and purchase histories. It
allows variance in the number of choice tasks per subject, and variance
in the number of choice alternatives per task.

Estimating a simple volumetric demand model is easy. Use the
`vd_est_vdm` function, and use at least 100,000 draws:

``` r
est_icecream <- icecream %>% vd_est_vdm(R=10000)
#> Using 16 cores
#>  MCMC in progress 
#> MCMC complete
#>  Total Time Elapsed: 0.15 minutes
```

Upper-level estimates can be summarized using `ec_estimates_MU`:

``` r
est_icecream %>% ec_estimates_MU()
#> # A tibble: 21 × 12
#>    attrib…¹ lvl   par      mean     sd `CI-5%` CI-95…² sig   model error refer…³
#>    <chr>    <chr> <chr>   <dbl>  <dbl>   <dbl>   <dbl> <lgl> <chr> <chr> <chr>  
#>  1 <NA>     <NA>  int   -3.22   0.545   -3.50  -2.71   TRUE  VD-c… EV1   <NA>   
#>  2 Brand    Blue… Bran… -0.695  0.155   -0.882 -0.458  TRUE  VD-c… EV1   BenNJe…
#>  3 Brand    Blue… Bran… -0.618  0.145   -0.796 -0.355  TRUE  VD-c… EV1   BenNJe…
#>  4 Brand    Brey… Bran… -0.0750 0.0894  -0.235  0.0546 FALSE VD-c… EV1   BenNJe…
#>  5 Brand    Drye… Bran… -0.546  0.126   -0.704 -0.342  TRUE  VD-c… EV1   BenNJe…
#>  6 Brand    Haag… Bran… -0.285  0.0861  -0.425 -0.150  TRUE  VD-c… EV1   BenNJe…
#>  7 Brand    Store Bran… -0.498  0.114   -0.650 -0.337  TRUE  VD-c… EV1   BenNJe…
#>  8 Flavor   Choc… Flav… -0.317  0.127   -0.493 -0.0434 TRUE  VD-c… EV1   Chocol…
#>  9 Flavor   Choc… Flav… -0.392  0.115   -0.563 -0.206  TRUE  VD-c… EV1   Chocol…
#> 10 Flavor   Cook… Flav… -0.415  0.111   -0.580 -0.235  TRUE  VD-c… EV1   Chocol…
#> # … with 11 more rows, 1 more variable: parameter <chr>, and abbreviated
#> #   variable names ¹​attribute, ²​`CI-95%`, ³​reference_lvl
```

Corresponding demand predictions can be obtained using the `vd_dem_vdm`
function. Here, we generate in-sample predictions:

``` r
dempres_icecream <-
  icecream %>%
  vd_dem_vdm(est = est_icecream)
#> Using 16 cores
#>  Computation in progress
```

The resulting output makes it easy to work with demand predictions
without obtaining posterior means too early. Demand prediction draws are
stored in a single column `.demdraws`.

``` r
dempres_icecream
#> # A tibble: 39,600 × 9
#>       id  task   alt     x     p Brand     Flavor      Size  .demdraws    
#>  * <int> <int> <int> <dbl> <dbl> <fct>     <fct>       <ord> <list>       
#>  1     1     1     1     8 0.998 Store     Neapolitan  16    <dbl [1,000]>
#>  2     1     1     2     0 0.748 Store     VanillaBean 16    <dbl [1,000]>
#>  3     1     1     3     0 1.25  BenNJerry Oreo        16    <dbl [1,000]>
#>  4     1     1     4     0 0.748 BenNJerry Neapolitan  16    <dbl [1,000]>
#>  5     1     1     5     0 2.49  HaagenDa  RockyRoad   4     <dbl [1,000]>
#>  6     1     1     6     0 1.25  HaagenDa  Oreo        16    <dbl [1,000]>
#>  7     1     1     7     0 1.12  BlueBunny Oreo        16    <dbl [1,000]>
#>  8     1     1     8     0 1.99  BlueBunny Neapolitan  4     <dbl [1,000]>
#>  9     1     1     9     0 0.622 Breyers   RockyRoad   16    <dbl [1,000]>
#> 10     1     1    10     0 3.49  Breyers   Vanilla     4     <dbl [1,000]>
#> # … with 39,590 more rows
```

We can aggregate (e.g., by subject `id`) using `ec_dem_aggregate`:

``` r
dempres_icecream %>% 
  ec_dem_aggregate('id')
#> # A tibble: 300 × 2
#>       id .demdraws    
#>    <int> <list>       
#>  1     1 <dbl [1,000]>
#>  2     2 <dbl [1,000]>
#>  3     3 <dbl [1,000]>
#>  4     4 <dbl [1,000]>
#>  5     5 <dbl [1,000]>
#>  6     6 <dbl [1,000]>
#>  7     7 <dbl [1,000]>
#>  8     8 <dbl [1,000]>
#>  9     9 <dbl [1,000]>
#> 10    10 <dbl [1,000]>
#> # … with 290 more rows
```

Once we have the desired aggregation level, we can obtain summaries
(e.g., posterior means) using `ec_dem_summarise`

``` r
dempres_icecream %>% 
  ec_dem_aggregate('id') %>%
  ec_dem_summarise()
#> # A tibble: 300 × 6
#>       id .demdraws     `E(demand)` `S(demand)` `CI-5%` `CI-95%`
#>    <int> <list>              <dbl>       <dbl>   <dbl>    <dbl>
#>  1     1 <dbl [1,000]>        39.8       13.5    19.8      64.6
#>  2     2 <dbl [1,000]>        99.6       27.0    57.8     145. 
#>  3     3 <dbl [1,000]>        30.3        5.79   20.8      40.1
#>  4     4 <dbl [1,000]>        89.7       27.4    50.4     139. 
#>  5     5 <dbl [1,000]>        31.8       17.5     9.57     65.8
#>  6     6 <dbl [1,000]>        15.9        8.76    4.64     31.6
#>  7     7 <dbl [1,000]>        73.2       20.0    49.9     108. 
#>  8     8 <dbl [1,000]>        51.2       20.3    23.9      91.6
#>  9     9 <dbl [1,000]>        13.9        4.43    7.05     21.7
#> 10    10 <dbl [1,000]>        37.5       11.1    19.2      55.8
#> # … with 290 more rows
```

Both `ec_dem_aggregate` and `ec_dem_summarise` simply apply common
`dplyr` and `purrr` functions.
