% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{get_metrics}
\alias{get_metrics}
\title{Get Model Metrics}
\usage{
get_metrics(
  mod,
  X = NULL,
  y = NULL,
  weights = NULL,
  threshold = NULL,
  threshold_type = c("response", "score")
)
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

\item{threshold}{Numeric vector of classification threshold values used to
calculate the accuracy, sensitivity, and specificity of the model. Defaults
to a range of risk probability thresholds from 0.1 to 0.9 by 0.1.}

\item{threshold_type}{Defines whether the \code{threshold} vector contains
risk probability values ("response") or threshold values expressed as scores
from the risk score model ("score"). Default: "response".}
}
\value{
Data frame with accuracy, sensitivity, and specificity for each threshold.
}
\description{
Calculates a risk model's accuracy, sensitivity, and specificity
given a set of data.
}
\examples{
y <- breastcancer[[1]]
X <- as.matrix(breastcancer[,2:ncol(breastcancer)])

mod <- risk_mod(X, y)
get_metrics(mod, X, y)

get_metrics(mod, X, y, threshold = c(150, 175, 200), threshold_type = "score")
}
