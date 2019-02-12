\name{srm.modelcompare}
\alias{srm.modelcompare}
\title{Comparing Models in Social Relations Model}
\usage{
srm.modelcompare(model.1, model.2, null = NULL, df = NULL, digits = 3)
}
\description{
Comparing Models in Social Relations Model.
}
\examples{
require(mlsrm)
data(srm_data)

# Null Model
m0 <- mlsrm(dv~1, "gid", "aid", "pid", srm_data)
# Prediction Model.1
m1 <- mlsrm(dv~iv, "gid", "aid", "pid", srm_data)
# Prediction Model.2
m2 <- mlsrm(dv~iv+m, "gid", "aid", "pid", srm_data)

# Comparing Model 1 with Model 2
srm.modelcompare(m1,m2,m0,df=1,digits=3)
}
