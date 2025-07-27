\name{GenerateCrossesfromElites}
\alias{GenerateCrossesfromElites}
\title{Generate crosses from elites}
\description{
Given a list of elite training sets, list of candidates the function makes npop new solutions by using crossover and mutation operators.
}
\usage{
GenerateCrossesfromElites(Elites, Candidates, npop, mutprob, mc.cores = 1,
                 mutintensity = 1, memoryfortabu = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Elites}{
a list of elite training sets
}
  \item{Candidates}{
a vector of identifiers of the individuals in the candidate set.
}
  \item{npop}{
number of training sets to generate.
}
  \item{mutprob}{
point mutation probability for each individual generated. Only one mutation per solution is allowed.  
}
 \item{mc.cores}{
number of cores to use.
}
 \item{mutintensity}{
mean of the poisson variable that is used to decide the number of mutations for each cross. 
}
 \item{memoryfortabu}{
tabu memory}
}

\value{
A list of npop training sets.
}

\author{
Deniz Akdemir
}
