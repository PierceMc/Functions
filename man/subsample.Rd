% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsamplepolygons.R
\name{subsample}
\alias{subsample}
\title{Return a rectangular sub-map of a polygon layer}
\usage{
subsample(sizex, sizey, polygons, rasterize = F, layer = "")
}
\arguments{
\item{sizex}{X axis length}

\item{sizey}{Y axis length}

\item{polygons}{Polygon layer to subsample}

\item{rasterize}{\link{optional} Return a raster object insteard of a polygon}

\item{layer}{\link{for rasterize} Layer of the polygon to convert to raster}
}
\value{
polygon or raster layer of size \link{siex, sizey} with the data from the input layer
}
\description{
Return a rectangular sub-map of a polygon layer
}
