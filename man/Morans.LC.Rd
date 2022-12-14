% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CostOfTravelMoransI.R
\name{Morans.LC}
\alias{Morans.LC}
\title{Calculate Morans I with a least cost path weighted matrix}
\usage{
Morans.LC(
  BaseDir,
  Variable,
  Values,
  dclasses,
  Species = "SBW",
  NSims = 32,
  Samples = 32,
  N = 2000,
  Transition = "transitionfunc",
  ncores = 4,
  verbose = F,
  SampleLength = 15
)
}
\arguments{
\item{Variable}{Name of the Variable in the filename}

\item{Values}{Vector of parameter values for the Moran's I to be calculated for}

\item{dclasses}{Distance classes for the correlogram}

\item{Species}{Species to use. Default=SBW}

\item{NSims}{How many model simulations. Default=32}

\item{Samples}{How many samples of the landscape to take. Default=32}

\item{N}{Points from the landscape to take each sample. Default=2000}

\item{Transition}{A function describing the cost of travel over non habitat patches. Default="transitionfunc"}

\item{ncores}{The number of cores available for calculation}

\item{verbose}{Boolean. Print status.}

\item{Basedir}{Directory containing the data}
}
\value{
List of 2. Data frame with Moran's I correlogram data and List of example landscapes.
}
\description{
Calculate Morans I with a least cost path weighted matrix
}
