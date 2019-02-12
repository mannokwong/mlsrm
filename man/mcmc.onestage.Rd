\name{mcmc.onestage}
\alias{mcmc.onestage}
\title{MCMC Simulations for One-Stage Model}
\usage{
mcmc.onestage(fittedModel, rep=20000)
}
\description{
MCMC Simulations for One-Stage Model.
This package is useful for testing
(1) simple slopes
(2) polynomial regression
}
\examples{
require(mlsrm)

# Simple slope tests
m1 <- mlsrm(M~IV+W+IVW, "gid", "aid", "pid", data)
mcmc.1 <- mcmc.onestage(m1)
mcmc.ci((mcmc.1[,2]+mcmc.1[,4])) # When W is high
mcmc.ci((mcmc.1[,2]-mcmc.1[,4])) # When W is low
}
