% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/echoice2.R
\name{vd_est_vdm_ssq}
\alias{vd_est_vdm_ssq}
\title{Estimate volumetric demand model accounting for set size variation (2nd order)}
\usage{
vd_est_vdm_ssq(vd, R=100000, keep=10)
}
\arguments{
\item{vd}{volumetric demand data (long format)}

\item{R}{draws}

\item{keep}{thinning}

\item{cores}{no of CPU cores to use (default: auto-detect)}

\item{control}{list containing additional settings}
}
\value{
est ec-draw object (List)
}
\description{
For more details on the model: https://dx.doi.org/10.2139/ssrn.3418383
This model REQUIRES variation in choice-set size
}
