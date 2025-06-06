% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic.R
\name{plot.cv_risk_mod}
\alias{plot.cv_risk_mod}
\title{Plot Risk Score Cross-Validation Results}
\usage{
\method{plot}{cv_risk_mod}(x, ...)
}
\arguments{
\item{x}{An object of class "cv_risk_mod", usually a result of a call to
\code{\link[=cv_risk_mod]{cv_risk_mod()}}.}

\item{...}{Additional arguments affecting the plot produced}
}
\value{
Object of class "ggplot".
}
\description{
Plots the mean auc for each \eqn{lambda_0} tested during cross-validation.
}
