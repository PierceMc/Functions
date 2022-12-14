% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.xyz.R
\name{as.xyz}
\alias{as.xyz}
\title{Converts a matrix of spatial values}
\usage{
as.xyz(d, raster = T)
}
\arguments{
\item{d}{Matrix of spatial values}

\item{raster}{Boolean. Inverts y values to account for raster standards}
}
\value{
Data frame of x position, y position, and values (z).
}
\description{
Takes a Matrix of spatial values and returns a
data frame of x, y, and z (values).
}
