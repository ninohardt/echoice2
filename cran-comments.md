## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
  -  This package uses OpenMP to speed up MCMC samplers. Please make sure it is part of your toolchain. This should not be an issue on most linux distributions and Windows with RTools installed. All examples are set to use 2 cores, but if set to NULL, all available cores will be used.