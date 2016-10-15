\name{SimReg3-package}
\alias{SimReg3-package}
\alias{SimReg3}
\docType{package}
\title{
Similarity Regression Functions
}
\description{
	Functions for performing Bayesian similarity regression,
    and evaluating the probability of association between sets of ontological terms
    and binary response vector. A random model is compared with one in which
    the log odds of a true response is linked to the semantic similarity
    between terms and a latent characteristic ontological profile.
}
\details{
\tabular{ll}{
Package: \tab SimReg3\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2016-10-14\cr
License: \tab GPL (>= 2)\cr
}
Key functions include sim_reg, for similarity regression of binary response variable against an ontologically encoded predictor. An example application would be inferring the probability of association between the presence of a rare genetic variant conditional on an ontologically encoded phenotype.
}
\author{
Daniel Greene <dg333@cam.ac.uk>

Maintainer: Daniel Greene <dg333@cam.ac.uk>
}
\references{
	D. Greene, NIHR BioResource, S. Richardson, E. Turro, `Phenotype similarity regression for identifying the genetic determinants of rare diseases', The American Journal of Human Genetics 98, 1-10, March 3, 2016.
}
\keyword{ HPO , MCMC, ontology }
