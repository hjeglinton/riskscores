% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic.R
\name{plot.risk_mod}
\alias{plot.risk_mod}
\title{Plot Risk Score Model Curve}
\usage{
\method{plot}{risk_mod}(x, score_min = NULL, score_max = NULL, ...)
}
\arguments{
\item{x}{An object of class "risk_mod", usually a result of a call to
\code{\link[=risk_mod]{risk_mod()}}.}

\item{score_min}{The minimum score displayed on the x-axis. The default is the
minimum score predicted from model's training data.}

\item{score_max}{The maximum score displayed on the x-axis. The default is the
maximum score predicted from model's training data.}

\item{...}{Additional arguments affecting the plot produced}
}
\value{
Object of class "ggplot".
}
\description{
Plots the linear regression equation associated with the integer risk score
model. Plots the scores on the x-axis and risk on the y-axis.
}
\examples{
y <- breastcancer[[1]]
X <- as.matrix(breastcancer[,2:ncol(breastcancer)])
mod <- risk_mod(X, y, lambda0 = 0.01)

plot(mod)
}
