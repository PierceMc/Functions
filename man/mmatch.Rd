% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmatch.R
\name{mmatch}
\alias{mmatch}
\title{Matches multiple columns in a dictionary to return a value}
\usage{
mmatch(from1, from2, to, datafrom1, datafrom2)
}
\arguments{
\item{from1}{Vector of from data}

\item{from2}{Vector of from data}

\item{to}{What the data will become}

\item{datafrom1}{vector containing values from from1 from the dataframe}

\item{datafrom2}{vector containing values from from2 from the dataframe}
}
\value{
Vector of target values based on from values
}
\description{
Like base match but accepts multiple 'from' arguments
}
