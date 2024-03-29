% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/echoice2.R
\name{dd_dem_sr}
\alias{dd_dem_sr}
\title{Discrete Choice Predictions (HMNL with attribute-based screening)}
\usage{
dd_dem_sr(dd, est, prob = FALSE, cores = NULL)
}
\arguments{
\item{dd}{data}

\item{est}{est}

\item{prob}{logical, report probabilities instead of demand}

\item{cores}{cores}
}
\value{
Draws of expected choice
}
\description{
Discrete Choice Predictions (HMNL with attribute-based screening)
}
\examples{
\donttest{
data(icecream_discrete)
icecream_est <- icecream_discrete \%>\% filter(id<20) \%>\% dd_est_hmnl_screen(R=10, cores=2)
#demand prediction
icecream_dempred <- icecream_discrete \%>\% filter(id<20) \%>\% 
 dd_dem_sr(icecream_est, cores=2)
 }

}
\seealso{
\code{\link[=dd_est_hmnl_screen]{dd_est_hmnl_screen()}} to generate demand predictions based on this model
}
