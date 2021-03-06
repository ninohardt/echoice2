% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/echoice2.R
\name{vd_est_vdm_screen}
\alias{vd_est_vdm_screen}
\title{Estimate volumetric demand model with attribute-based conjunctive screening}
\usage{
vd_est_vdm_screen(vd, R=100000, keep=10)
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
See https://dx.doi.org/10.2139/ssrn.2770025 for more details
}
\examples{
\dontrun{
icecream_est <- icecream \%>\% vd_est_vdm_screen(R=50000)
}
}
