% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dGAselID.R
\name{RandomizePop}
\alias{RandomizePop}
\title{RandomizePop}
\usage{
RandomizePop(population)
}
\arguments{
\item{population}{Population of chromosome sets in current generation.}
}
\description{
Generates a random population for the next generation.
}
\examples{
\dontrun{
library(ALL)
data(ALL)

demoALL<-ALL[1:12,1:8]

population01<-InitialPopulation(demoALL, 4, 4)
population01
RandomizePop(population01)
 }
}

