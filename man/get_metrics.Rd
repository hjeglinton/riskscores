% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{get_metrics}
\alias{get_metrics}
\title{Get Model Metrics}
\usage{
get_metrics(mod, X = NULL, y = NULL, weights = NULL)
}
\arguments{
\item{mod}{An object of class \code{risk_mod}, usually a result of a call to
\code{\link[=risk_mod]{risk_mod()}}.}

\item{X}{Input covariate matrix with dimension \eqn{n \times p};
every row is an observation.}

\item{y}{Numeric vector for the (binomial) response variable.}

\item{weights}{Numeric vector of length \eqn{n} with weights for each
observation. Unless otherwise specified, default will give equal weight to
each observation.}
}
\value{
List with deviance (dev), accuracy (acc), sensitivity (sens), and
specificity (spec).
}
\description{
Calculates a risk model's deviance, accuracy, sensitivity, and specificity
given a set of data.
}
\examples{
y <- breastcancer[[1]]
X <- as.matrix(breastcancer[,2:ncol(breastcancer)])

mod <- risk_mod(X, y)
get_metrics(mod, X, y)
}
