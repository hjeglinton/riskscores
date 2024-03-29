% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{stratify_folds}
\alias{stratify_folds}
\title{Generate Stratified Fold IDs}
\usage{
stratify_folds(y, nfolds = 10, seed = NULL)
}
\arguments{
\item{y}{Numeric vector for the (binomial) response variable.}

\item{nfolds}{Number of folds (default: 10).}

\item{seed}{An integer that is used as argument by \code{set.seed()} for
offsetting the random number generator. Default is to not set a
particular randomization seed.}
}
\value{
Numeric vector with the same length as \code{y}.
}
\description{
Returns a vector of fold IDs that preserves class proportions.
}
\examples{
y <- rbinom(100, 1, 0.3)
foldids <- stratify_folds(y, nfolds = 5)
table(y, foldids)
}
