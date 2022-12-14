% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Empirical.R
\name{CoxDataPrep}
\alias{CoxDataPrep}
\title{Prepare output of main.py for cox analysis}
\usage{
CoxDataPrep(
  x,
  finalyear = 2020,
  strata = 23000,
  OnsetDefol = 50,
  Mortdefol = 600
)
}
\arguments{
\item{x}{Data from main.py}

\item{finalyear}{Final year of data collection. Default=2020}

\item{strata}{Size of bins for stratification in metres. Default=23000}
}
\value{
Year of first defoliation
}
\description{
Prepare output of main.py for cox analysis
}
