% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/echoice2.R
\name{ec_estimates_SIGMA_corr}
\alias{ec_estimates_SIGMA_corr}
\title{Obtain posterior mean estimates of upper level correlations}
\usage{
ec_estimates_SIGMA_corr(est)
}
\arguments{
\item{est}{is an echoice draw object (list)}
}
\value{
estimates of upper level correlations
}
\description{
Obtain posterior mean estimates of upper level correlations
}
\examples{
\dontrun{
#obtain correlations and visualize using corrr package
est \%>\% ec_estimates_SIGMA_corr \%>\% 
corrr::as_cordf() \%>\% 
corrr::rplot(print_cor = TRUE, legend = FALSE)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
}
