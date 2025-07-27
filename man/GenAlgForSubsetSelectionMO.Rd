\name{GenAlgForSubsetSelectionMO}
\alias{GenAlgForSubsetSelectionMO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Genetic algorithm for subset selection no given test with multiple criteria for Multi Objective Optimized Experimantal Design.}
\description{
It uses a nondominated selection genetic algorithm to find the solutions on the frontier that optimizes several design criteria at the same time. The test set is taken as the complement of the training individuals. 
}

\usage{
GenAlgForSubsetSelectionMO(Pcs = NULL, Dist = NULL, Kernel = NULL, Candidates, Test,
                 ntoselect, selectionstats, selectionstatstypes,
                 plotdirections, npopGA = 100, mutprob = 0.8,
                 mutintensity = 1, nitGA = 500, lambda = 1e-06,
                 plotiters = FALSE, mc.cores = 1, InitPop = NULL, C =
                 NULL, axes.labels = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pcs}{
Principal components matrix for the individuals
}
 \item{Dist}{
Distance matrix for the individuals
}
 \item{Kernel}{
Kernel matrix for the individuals
}

 \item{Candidates}{The set of individuals from which the training set is selected.
}
\item{Test}{The set of individuals for which the predictions based on the model built using the training set are needed. 
}
   \item{ntoselect}{number of individuals to select in the training set.
}

 \item{selectionstats}{a vector of design optimization criteria
}

 \item{selectionstatstypes}{a vector describing the type of optimality criteria used in \code{selectionstats}.
}

\item{plotdirections}{A vector that is used to change the sign of statistics while plotting. This doesnt effect the optimization.
}
 
  \item{npopGA}{
 genetic algorithm parameter, number of solutions at each iteration
}
  \item{mutprob}{
 genetic algorithm parameter, probability of mutation for each generated solution.
}
 \item{mutintensity}{
genetic algorithm parameter, mean of the poisson variable that is used to decide the number of mutations for each cross. 
}
  \item{nitGA}{
 genetic algorithm parameter, number of iterations.  
}
\item{lambda}{
 scalar shrinkage parameter (\eqn{\lambda>0}).
}
\item{plotiters}{
 plot the convergence: TRUE or FALSE. Default is TRUE.
}

\item{mc.cores}{
 number of cores to use.
}
\item{InitPop}{
a list of initial solutions
}

 \item{C}{
Contrast Matrix.
}

\item{axes.labels}{Labels for the axes for plotting iterations}
}

\value{
A list of length 2. The first item in the list is the list of solutions found by the algorithm. The second item is a matrix of criteria values for the solutions in the first list. 
}

\author{
Deniz Akdemir
}


\examples{
\dontrun{
library(STPGA)
library(GenomicMating)

data(WheatData)


Msvd<-svd(scale(Wheat.M, scale=F, center=T), nu=50, nv=50)
Dgeno<-as.matrix(dist(scale(Wheat.M, scale=F, center=T)))^2
P<-Wheat.M\%*\%Msvd$v
dim(P)
rownames(Dgeno)<-colnames(Dgeno)<-rownames(P)<-rownames(Wheat.M)
test<-sample(rownames(P), 25)
candidates<-setdiff(rownames(P), test)
outnewprog<-GenAlgForSubsetSelectionMO(Pcs=P,Dist=Dgeno,
Candidates=candidates,Test=test,ntoselect=75,
selectionstats=list("DOPT", "neg_dist_in_train2", "dist_to_test2"),
selectionstatstypes=c("Pcs", "Dist", "Dist"),
plotdirections=c(1,1,1),
npopGA=300,  mutprob=1, mutintensity=2,
nitGA=100, plotiters=TRUE, mc.cores=1, InitPop=NULL)

#####Best solution according to ideal solution concept
outnewprog[[1]][[which.min(disttoideal(outnewprog[[2]]))]]
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

