\name{legacy-compatibility}
\alias{legacy-compatibility}
\alias{Amat.pieces}
\alias{calculate_population_distances}
\alias{dist_to_test}
\alias{dist_to_test2}
\alias{disttoideal}
\alias{evaluate_population_legacy}
\alias{evaluate_population_optimized}
\alias{GenAlgForSubsetSelectionMO}
\alias{GenAlgForSubsetSelectionMONoTest}
\alias{neg_dist_in_train}
\alias{neg_dist_in_train2}
\alias{unified_distance_criterion}
\title{Legacy compatibility wrappers}
\usage{
Amat.pieces(M, pieces = 10, mc.cores = 1)
calculate_population_distances(solutions, method = c("hamming", "jaccard", "sorensen"))
dist_to_test(train, test, Dst, lambda = NULL, C = NULL)
dist_to_test2(train, test, Dst, lambda = NULL, C = NULL)
disttoideal(X, method = "euclidean", handle_zeros = "warning")
evaluate_population_legacy(
  population,
  errorstat,
  Test,
  P,
  K = NULL,
  lambda = 1e-06,
  C = NULL,
  mc.cores = 1,
  Vg = NULL,
  Ve = NULL
)
evaluate_population_optimized(
  population,
  P,
  Test,
  errorstat = "PEVMEAN2",
  lambda = 1e-06,
  C = NULL,
  mc.cores = 1,
  K = NULL,
  Vg = NULL,
  Ve = NULL
)
GenAlgForSubsetSelectionMO(...)
GenAlgForSubsetSelectionMONoTest(...)
neg_dist_in_train(train, test = NULL, Dst, lambda = NULL, C = NULL)
neg_dist_in_train2(train, test = NULL, Dst, lambda = NULL, C = NULL)
unified_distance_criterion(
  train,
  test,
  Dst,
  criterion = c("max_to_test", "mean_to_test", "neg_min_internal", "neg_mean_internal"),
  lambda = NULL,
  C = NULL
)
}
\description{
Deprecated or compatibility wrappers retained for older STPGA code.
New code should use the corresponding modern functions documented elsewhere,
including \code{compute_amatrix()}, \code{population_distances()},
\code{distance_train_to_test_max()}, \code{distance_train_to_test_mean()},
\code{distance_internal_min()}, \code{distance_internal_mean()},
\code{distance_criterion()}, \code{evaluate_population()},
\code{subset_ga_multiobjective()}, and
\code{subset_ga_multiobjective_single()}.
}
\arguments{
\item{M, P, X}{Input matrix.}
\item{pieces}{Number of marker chunks for relationship matrix computation.}
\item{mc.cores}{Number of parallel worker cores.}
\item{solutions, population}{Candidate solutions to evaluate.}
\item{method}{Distance method.}
\item{train, test, Test}{Training and test set identifiers.}
\item{Dst}{Distance matrix.}
\item{lambda}{Ridge parameter.}
\item{C}{Optional contrast matrix.}
\item{handle_zeros}{How zero ranges are handled.}
\item{errorstat}{Legacy criterion name.}
\item{K}{Optional kinship matrix.}
\item{Vg, Ve}{Genetic and residual variance components.}
\item{criterion}{Distance criterion name.}
\item{...}{Arguments passed to the corresponding modern implementation.}
}
\value{
The value returned by the corresponding modern implementation.
}
