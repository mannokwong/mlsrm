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
m1 <- mlsrm(m~iv, "gid", "aid", "pid", srm_data)
m2 <- mlsrm(dv~iv+m, "gid", "aid", "pid", srm_data)
mcmc.1 <- mcmc.twostage(m1, m2)
mcmc.ci(mcmc.1[,2]*mcmc.1[,5])

}
