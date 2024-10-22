% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpit.R
\name{cpdf}
\alias{cpdf}
\title{Conditional PDF}
\usage{
cpdf(object, newdata, cores = 1)
}
\arguments{
\item{object}{an object of class \code{cvinereg}.}

\item{newdata}{matrix of response and covariate values for which to compute
the conditional density}

\item{cores}{integer; the number of cores to use for computations.}
}
\description{
Calculates the conditional density of the response given the covariates.
}
\examples{
# simulate data
x <- matrix(rnorm(200), 100, 2)
y <- x \%*\% c(1, -2)
data <- data.frame(y = y, x = x)

# fit vine regression model
fit <- cvinereg(y ~ ., data)
# pdf
d <- cpdf(fit, data)
# cll
sum(log(d))

}
