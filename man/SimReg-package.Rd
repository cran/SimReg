\name{SimReg-package}
\alias{SimReg-package}
\alias{SimReg}
\docType{package}
\title{
Similarity Regression Functions
}
\description{
Functions for performing Bayesian `similarity regression',
specialised for modelling the association between a HPO-encoded phenotype
and binary genotype. A `no association' model is compared with one in which
the log odds of having a rare genotype is linked to the semantic similarity
between patient phenotype and a latent characteristic phenotype. The method
estimates the probability of an association together with an HPO-coded phenotype
characteristic of the disease. 
}
\details{
\tabular{ll}{
Package: \tab SimReg\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2016-01-12\cr
License: \tab GPL (>= 2)\cr
}
Key functions include sim.reg, for performing similarity regression of binary response variable against an ontologically encoded predictor. The typical use would be for fitting a model for the presence of a rare genetic variant conditional on an ontologically encoded phenotype.
}
\author{
Daniel Greene <dg333@cam.ac.uk>

Maintainer: Daniel Greene <dg333@cam.ac.uk>
}
\references{
	D. Greene, NIHR BioResource, S. Richardson, E. Turro, `Phenotype similarity regression for identifying the genetic determinants of rare diseases', The American Journal of Human Genetics 98, 1-10, March 3, 2016.
}
\keyword{ HPO , MCMC, ontology }
