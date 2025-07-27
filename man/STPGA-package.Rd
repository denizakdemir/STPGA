\name{STPGA-package}
\alias{STPGA-package}
\alias{STPGA}
\docType{package}
\title{Selection of Training Populations by Genetic Algorithm}
\description{This package can be utilized to select a (test data) calibrated training population in high dimensional prediction problems. More specifically, the package contains a genetic algorithm that tries to minimize a design criterion defined for subsets of a certain size selected from a larger set. 
 }
\details{
\tabular{ll}{
Package: \tab STPGA\cr
Type: \tab Package\cr
Version: \tab 5.0\cr
Date: \tab 2018-07-20\cr
License: \tab GPL-3\cr
}

The package is  useful for high dimensional prediction problems where per individual cost of observing / analyzing the response variable is high and therefore a small number of training examples is sought or when the candidate set from which the training set must be chosen (is not representative of the test data set). 

The function "GenAlgForSubsetSelection" uses a simple genetic algorithm to identify a training set of a specified size from a larger set of candidates which minimizes an optimization criterion for a known test set. The function "GenAlgForSubsetSelectionNoTest" tries to identify a training set of a specified size from a larger set of candidates which minimizes an optimization criterion.
}
\author{
Maintainer: Deniz Akdemir <deniz.akdemir.work@gmail.com>
}
\references{
References: 
Akdemir, Deniz. "Training population selection for (breeding value) prediction." arXiv preprint arXiv:1401.7953 (2014).
}
