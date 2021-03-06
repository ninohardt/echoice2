% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/echoice2.R
\name{ec_estimates_screen}
\alias{ec_estimates_screen}
\title{Summarize attribute-based screening parameters}
\usage{
ec_estimates_screen(est, quantiles=c(.05,.95))
}
\arguments{
\item{est}{is an echoice draw object (list)}

\item{quantiles}{quantile for CI}
}
\value{
tibble with screening summaries
}
\description{
Summarize attribute-based screening parameters
}
