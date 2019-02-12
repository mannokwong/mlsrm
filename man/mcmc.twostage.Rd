\name{mcmc.twostage}
\alias{mcmc.twostage}
\title{MCMC Simulations for Two-Stage Model}
\usage{
mcmc.twostage(stage1_model,stage2_model, rep=20000)
}
\description{
MCMC Simulations for Two-Stage Model.
This package is useful for testing
(1) indirect effect
(2) conditional indirect effect
(3) polynomial regression (with moderated mediation)
}
\examples{
require(mlsrm)
...

# Testing indirect effect
m1 <- mlsrm(M~IV, "GID", "AID", "PID", data)
m2 <- mlsrm(DV~IV+M, "GID", "AID", "PID", data)
mcmc.1 <- mcmc.twostage(m1, m2)
mcmc.ci(mcmc.1[,2]*mcmc.1[,5])

# Testing conditional indirect effect
m1 <- mlsrm(M~IV+W+IVW, "GID", "AID", "PID", data)
m2 <- mlsrm(DV~IV+W+IVW+M, "GID", "AID", "PID", data)
mcmc.2 <- mcmc.twostage(m1, m2)
mcmc.ci((mcmc.2[,2]+mcmc.2[,4])*mcmc.2[,9]) # When W is high
mcmc.ci((mcmc.2[,2]-mcmc.2[,4])*mcmc.2[,9]) # When W is low


}
