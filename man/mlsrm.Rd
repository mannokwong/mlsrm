\name{srm}
\alias{srm}
\title{Social Relations Model}
\usage{
mlsrm(formula, gid, aid, pid, data)
}
\description{
Social Relations Model.
}
\examples{
require(mlsrm)

# Reading SPSS files (.sav)
require(foreign)
# for Window users
data = data.frame(read.spss(choose.files()))
# for Mac users
data = data.frame(read.spss(file.choose()))

# Null Model
m0 <- mlsrm(DV~1, "GID", "AID", "PID", data)
# Prediction Model.1
m1 <- mlsrm(DV~IV1, "GID", "AID", "PID", data)
# Prediction Model.2
m2 <- mlsrm(DV~IV1+IV2, "GID", "AID", "PID", data)
}
