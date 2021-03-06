% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wu.R
\name{wu.fit}
\alias{wu.fit}
\title{Fit spike and slab Cox PHMs}
\usage{
wu.fit(Y, delta, X)
}
\arguments{
\item{Y}{survival times}

\item{delta}{censoring indicator, 0: uncensored, 1: censored}

\item{X}{design matrix}
}
\value{
list with coefficient values and number of non-zero coefficients
}
\description{
Fit spike and slab Cox PHMs
}
\examples{
# Generate data
set.seed(1)
n <- 250; omega <- 1; censoring_lvl <- 0.4

b <- c(1, 1, rep(0, 250))
p <- length(b)
X <- matrix(rnorm(n * p), nrow=n)
y <- runif(nrow(X))
Y <- log(1 - y) / - (exp(X \%*\% b) * omega)

delta  <- runif(n) > censoring_lvl   # 0: censored, 1: uncensored
Y[!delta] <- Y[!delta] * runif(sum(!delta))

fit <- wu.fit(Y, delta, X)

}
