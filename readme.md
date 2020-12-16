echoice2
================

# ![echoice2 0.1](man/figures/echoicelogo.png)

Package website: [ninohardt.de/echoice2](http://ninohardt.de/echoice2)

This package contains choice models with economic foundation.

For more theoretical background, please refer to the chapter ‘economic
foundations of conjoint’ in the follwing handbook:
<https://www.elsevier.com/books/handbook-of-the-economics-of-marketing/dube/978-0-444-63759-8>

For details on how to use echoice, please refer to the vignette in the
package. It illustrates model estimation, predictions and model
evaluation.

echoice2 is a complete rewrite of the original echoice package. All key
functions are written in c++ and use openMP for multi-threaded
computing. echoice2 follows tidy principles and integrated nicely with
dplyr. It can be used to generate choice volume/share simulators, though
no frot-end is built into the package yet.

### Functionality

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

### Updates

The original echoice package is still available here:
<https://bitbucket.org/ninohardt/echoice>

### Installation

This package has not been published on CRAN yet. You can install it
using the remotes package.

``` r
remotes::install_github("ninohardt/echoice2")
```

Please make sure you are using a recent version of R (4.0.x) and you
update packages such as tidyverse, Rcpp and RcppArmadillo.

If you are installing this on windows, make you install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) first.

If you are using OSX, make you have CLI, XQuartz and potentially other
things installed.

If you are using linux, this should work.

Binaries for windows can be downloaded
[here](http://ninohardt.de/echoice2/echoice2_0.1.zip "echoice2")

### Using echoice

Please read the vignette. It illustrates a complete workflow from
estimation to market simulation. The vignette can also be found on the
[package website](http://ninohardt.de/echoice2/articles/echoice2.html).
