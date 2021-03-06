% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/echoice2.R
\name{ec_screen_summarise}
\alias{ec_screen_summarise}
\title{Summarize posterior draws of screening}
\usage{
ec_dem_summarise(de,quantiles)
}
\arguments{
\item{sc}{tibble containing screening draws in .screendraws}

\item{quantiles}{Quantiles for Credibility Intervals (default: 90\% interval)}
}
\value{
Summary of screening probabilities
}
\description{
Adds summaries of posterior draws of demand to tibble.
(using the new demand draw format)
}
