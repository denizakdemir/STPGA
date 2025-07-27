\name{disttoideal}

\alias{disttoideal}



\title{Calculate the distance of solutions from the 'ideal' solution.}
\description{
This function calculates the distance of X to the vector of the minimums of columns of X after transforming the variables in X to the interval \eqn{[0,1]}.
}
\usage{
disttoideal(X)
}
\arguments{
    \item{X}{
A matrix of the criteria values. One solution each row, columns correspond to the different criteria.
}
}
\value{Vector of distances elements corresponding to each row of X}

\author{Deniz Akdemir}

