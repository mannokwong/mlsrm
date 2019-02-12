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

data(srm_data)

# Null Model
m0 <- mlsrm(dv~1, "gid", "aid", "pid", srm_data)
# Prediction Model.1
m1 <- mlsrm(dv~iv, "gid", "aid", "pid", srm_data)
# Prediction Model.2
m2 <- mlsrm(dv~iv+m, "gid", "aid", "pid", srm_data)
}
