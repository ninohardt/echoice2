
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
  [**0.2.5**](https://github.com/ninohardt/echoice2/releases/tag/v0.2.5):

  - removed stringr dependency
    
- Version
  [**0.2.4**](https://github.com/ninohardt/echoice2/releases/tag/v0.2.4):

  - A lot of work because of Roxygen’s [breaking
    change](https://github.com/r-lib/roxygen2/issues/1491), from
    recommended use of docType to breaking it
  - Nothing changed functionality-wise

- Version
  [**0.2.3**](https://github.com/ninohardt/echoice2/releases/tag/v0.2.3):

  - The vignette “Importing list-of-lists choice data and discrete
    choice modeling with echoice2” demonstrates how to import
    list-of-lists style choice data, convert it for use with echoice2
    and fit choice models. Related convenience functions have been added
    to the package.
  - CRAN release to follow soon

- Version
  [**0.2.2**](https://github.com/ninohardt/echoice2/releases/tag/v0.2.2):

  - Added a vignette on volumetric choice modeling with and without
    conjunctive screening

- Version
  [**0.2.1**](https://github.com/ninohardt/echoice2/releases/tag/v0.2.1):

  - Package compiles in absense of OpenMP Support
  - MacOS CRAN binaries are compiled without OpenMP - for full speed
    compile yourself or use binaries supplied here

- Version
  [**0.2.0**](https://github.com/ninohardt/echoice2/releases/tag/v0.20):

  - Initial CRAN release
  - No new functionality, but cleaner code

- Version
  [**0.1.6**](https://github.com/ninohardt/echoice2/releases/tag/v0.16):

  - Stability/Performance improvements for demand simulators

- Version
  [**0.1.5**](https://github.com/ninohardt/echoice2/releases/tag/v0.15):

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
#>  Total Time Elapsed: 0.17 minutes
```

Upper-level estimates can be summarized using `ec_estimates_MU`:

``` r
est_icecream %>% ec_estimates_MU()
#> # A tibble: 21 × 12
#>    attribute lvl         par      mean     sd `CI-5%` `CI-95%` sig   model error
#>    <chr>     <chr>       <chr>   <dbl>  <dbl>   <dbl>    <dbl> <lgl> <chr> <chr>
#>  1 <NA>      <NA>        int    -3.22  0.547   -3.56   -2.52   TRUE  VD-c… EV1  
#>  2 Brand     BlueBell    Brand… -0.713 0.169   -0.944  -0.471  TRUE  VD-c… EV1  
#>  3 Brand     BlueBunny   Brand… -0.731 0.169   -0.945  -0.394  TRUE  VD-c… EV1  
#>  4 Brand     Breyers     Brand… -0.117 0.0976  -0.295   0.0331 FALSE VD-c… EV1  
#>  5 Brand     Dryers      Brand… -0.554 0.122   -0.705  -0.358  TRUE  VD-c… EV1  
#>  6 Brand     HaagenDa    Brand… -0.358 0.0937  -0.510  -0.211  TRUE  VD-c… EV1  
#>  7 Brand     Store       Brand… -0.526 0.126   -0.707  -0.348  TRUE  VD-c… EV1  
#>  8 Flavor    ChocChip    Flavo… -0.393 0.113   -0.566  -0.213  TRUE  VD-c… EV1  
#>  9 Flavor    ChocDough   Flavo… -0.435 0.128   -0.618  -0.192  TRUE  VD-c… EV1  
#> 10 Flavor    CookieCream Flavo… -0.443 0.111   -0.611  -0.256  TRUE  VD-c… EV1  
#> # ℹ 11 more rows
#> # ℹ 2 more variables: reference_lvl <chr>, parameter <chr>
```

Corresponding demand predictions can be obtained using the `vd_dem_vdm`
function. Here, we generate in-sample predictions:

``` r
dempres_icecream <-
  icecream %>%
  vd_dem_vdm(est = est_icecream)
#> Using 16 cores
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
#> # ℹ 39,590 more rows
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
#> # ℹ 290 more rows
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
#>  1     1 <dbl [1,000]>        39.5       12.6    20.7      62.5
#>  2     2 <dbl [1,000]>        99.1       27.5    56.7     146. 
#>  3     3 <dbl [1,000]>        31.6        6.05   21.7      41.5
#>  4     4 <dbl [1,000]>        87.9       28.9    48.0     138. 
#>  5     5 <dbl [1,000]>        32.2       17.9    10.7      68.2
#>  6     6 <dbl [1,000]>        16.1        8.78    4.96     33.6
#>  7     7 <dbl [1,000]>        72.1       22.4    47.7     110. 
#>  8     8 <dbl [1,000]>        49.7       20.5    21.5      88.9
#>  9     9 <dbl [1,000]>        13.7        4.66    6.48     22.3
#> 10    10 <dbl [1,000]>        38.0       11.5    19.3      56.7
#> # ℹ 290 more rows
```

Both `ec_dem_aggregate` and `ec_dem_summarise` simply apply common
`dplyr` and `purrr` functions.
