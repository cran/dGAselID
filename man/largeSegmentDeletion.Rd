% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dGAselID.R
\name{largeSegmentDeletion}
\alias{largeSegmentDeletion}
\title{largeSegmentDeletion}
\usage{
largeSegmentDeletion(individuals, chrConf, mutationChance)
}
\arguments{
\item{individuals}{dataset returned by Individuals().}

\item{chrConf}{Configuration of chromosomes returned by splitChromosomes().}

\item{mutationChance}{Chance for a large segment deletion mutation to occur.}
}
\description{
Operator for the large segment deletion.
}
\examples{
\dontrun{
library(ALL)
data(ALL)

demoALL<-ALL[1:12,1:8]

set.seed(1234)
population<-InitialPopulation(demoALL, 4, 9)
individuals<-Individuals(population)

chrConf<-splitChromosomes(demoALL, 2)
chrConf
individuals

set.seed(123)
largeSegmentDeletion(individuals, chrConf, 20)
 }
}

