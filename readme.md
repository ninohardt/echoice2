
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

For details and examples on how to use echoice, please refer to the
[vignette](http://ninohardt.de/echoice2/articles/echoice2.html) in the
package. It illustrates model estimation, predictions and model
evaluation.

All key functions are written in c++ and use
[openMP](https://www.openmp.org) for multi-threaded computing. C++
integration in R is facilitated by [Rcpp](http://www.rcpp.org).

`echoice2` (largely) follows tidy principles and integrated nicely with
`dplyr`. It can be used to generate choice volume/share simulators,
though no front-end is built into the package yet.

`echoice2` is a complete rewrite of the original `echoice` package,
which is still available on
[Bitbucket](https://bitbucket.org/ninohardt/echoice).

## News

-   **Looking for co-developers!**

-   Next: upper-level covariates, effects-coding, discrete and
    continuous attributes, discrete choice example vignette, …

-   Version
    [**0.16**](https://github.com/ninohardt/echoice2/releases/tag/v0.16):

    -   Stability/Performance improvements for demand simulators

-   Version
    [**0.15**](https://github.com/ninohardt/echoice2/releases/tag/v0.15):

    -   Initial release
    -   faster and more efficient screening model estimation
    -   improved demand predictions: posterior demand draws are now
        stored in a single column - this is a *major* improvement for
        dealing with demand predictions!
    -   some bug-fixes and documentation improvements

## Installation

``` r
# install.packages("remotes")

#Turn off warning-error-conversion, because the tiniest warning stops installation
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

#install from github
remotes::install_github("ninohardt/echoice2")

#install from github and compile/include vingette
remotes::install_github("ninohardt/echoice2", build_vignettes = TRUE)
```

### Installation notes

-   Make sure you are using a recent R4.0.x version, `tidyerse`, `Rcpp`
    and `RcppArmadillo` before building the package from Github.

-   If you use Linux it should just work. If not, you might need to
    install libgomp before compilation.

-   If you are using OSX, you may have to install CLI, XQuartz and
    potentially other things that Apple removed from OSX. Google ‘+OSX
    +Rcpp’ if you run into trouble.

-   If you are using Windows, install
    \[Rtools\](<https://cran.r-project.org/bin/windows/Rtools/)> first.

### Binaries

Binaries are made available with
[releases](https://github.com/ninohardt/echoice2/releases) (under
‘Assets’).

You can install the binary version from the ‘Packages’ tab in RStudio.
Select ‘Install’ and change ‘Install from’ to ‘Package Archive File’. Or
you use `install.packages` on the command line and point it to the
downloaded .zip file.

## Functionality

The following models are implemented (including estimation and
prediction):

-   Discrete Choice (HMNL)
-   Discrete Choice, attribute-based screening (not including price)
-   Discrete Choice, attribute-based screening (including price)
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
-   attributes defining the choice alternatives (factor, and soon
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
#>  Total Time Elapsed: 0.15 minutes
```

Upper-level estimates can be summarized using `ec_estimates_MU`:

``` r
est_icecream %>% ec_estimates_MU()
#> # A tibble: 21 x 12
#>    attribute lvl   par      mean     sd `CI-5%` `CI-95%` sig   model error
#>    <chr>     <chr> <chr>   <dbl>  <dbl>   <dbl>    <dbl> <lgl> <chr> <chr>
#>  1 <NA>      <NA>  int   -3.24   0.550   -3.55   -2.60   TRUE  VD-c~ EV1  
#>  2 Brand     Blue~ Bran~ -0.678  0.162   -0.891  -0.419  TRUE  VD-c~ EV1  
#>  3 Brand     Blue~ Bran~ -0.605  0.151   -0.812  -0.405  TRUE  VD-c~ EV1  
#>  4 Brand     Brey~ Bran~ -0.0993 0.0963  -0.264   0.0584 FALSE VD-c~ EV1  
#>  5 Brand     Drye~ Bran~ -0.561  0.122   -0.721  -0.402  TRUE  VD-c~ EV1  
#>  6 Brand     Haag~ Bran~ -0.324  0.0880  -0.452  -0.188  TRUE  VD-c~ EV1  
#>  7 Brand     Store Bran~ -0.478  0.117   -0.646  -0.303  TRUE  VD-c~ EV1  
#>  8 Flavor    Choc~ Flav~ -0.393  0.132   -0.619  -0.195  TRUE  VD-c~ EV1  
#>  9 Flavor    Choc~ Flav~ -0.437  0.129   -0.643  -0.229  TRUE  VD-c~ EV1  
#> 10 Flavor    Cook~ Flav~ -0.393  0.101   -0.535  -0.227  TRUE  VD-c~ EV1  
#> # ... with 11 more rows, and 2 more variables: reference_lvl <chr>,
#> #   parameter <chr>
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
#> # A tibble: 39,600 x 9
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
#> # ... with 39,590 more rows
```

We can aggregate (e.g., by subject `id`) using `ec_dem_aggregate`:

``` r
dempres_icecream %>% 
  ec_dem_aggregate('id')
#> # A tibble: 300 x 2
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
#> # ... with 290 more rows
```

Once we have the desired aggregation level, we can obtain summaries
(e.g., posterior means) using `ec_dem_summarise`

``` r
dempres_icecream %>% 
  ec_dem_aggregate('id') %>%
  ec_dem_summarise()
#> # A tibble: 300 x 6
#>       id .demdraws     `E(demand)` `S(demand)` `CI-5%` `CI-95%`
#>    <int> <list>              <dbl>       <dbl>   <dbl>    <dbl>
#>  1     1 <dbl [1,000]>        40.2       14.1    19.7      66.7
#>  2     2 <dbl [1,000]>        99.3       27.7    55.5     148. 
#>  3     3 <dbl [1,000]>        31.0        6.04   21.9      41.7
#>  4     4 <dbl [1,000]>        85.9       28.3    44.8     137. 
#>  5     5 <dbl [1,000]>        32.7       17.8     9.93     67.9
#>  6     6 <dbl [1,000]>        15.5        8.23    4.65     30.0
#>  7     7 <dbl [1,000]>        72.9       22.2    48.5     105. 
#>  8     8 <dbl [1,000]>        52.8       19.9    25.2      89.7
#>  9     9 <dbl [1,000]>        13.8        4.75    6.46     22.1
#> 10    10 <dbl [1,000]>        38.1       10.8    20.2      55.8
#> # ... with 290 more rows
```

Both `ec_dem_aggregate` and `ec_dem_summarise` simply apply common
`dplyr` and `purrr` functions.
