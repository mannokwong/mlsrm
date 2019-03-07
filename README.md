# mlsrm
**mlsrm** enables the estimation of Social Relations Model using a multilevel approach described by Snijders &amp; Kenny (1999). The code (pdSRM.R) was originally developed by Andrew P Knight (knightap@wustl.edu) in 2016, and I modified and extended his code to make the analysis easier.

## Version History
+ **1.0** (7/3/2019)  
I have created an all-in-one function **mlsrm()** to make the social relations analysis easier.  
Other important features include:  
**mcmc.onestage()**, **mcmc.twostage()** and **mcmc.ci()** - Observed data is used to estimate the model parameters and their asymptotic covariance. These are then drawn to derive the CIs for each parameter of interest using Monte Carlo simulation. These functions are useful for assessing mediation and response surface analysis.  
**srm.modelcompare()** - Comparing fit between two models.

## Instructions
Install and Load Package
Open an R console or RStudio window. (R can be downloaded for free from https://cran.r-project.org; RStudio can be downloaded for free from https://www.rstudio.com/)

Install R package "mlsrm" through Github by pasting and running the following commands in R console or RStudio:
```R
install.packages("devtools") 
library("devtools") 
install_github("mannokwong/mlsrm")  
```

To use mlsrm, you need to call the package first:
```R
library("mlsrm") 
```

## Example Implementation
Prepare the dataset in long format (example is given below): 

ID|gid|aid|pid|iv|m|dv
-|-|-|-|-|-|-
1|10|11|12|1|5|5
2|10|11|13|1|5|5
3|10|11|14|1|5|5
4|10|11|15|1|5|5
5|10|12|11|2|4|3
...|...|...|...|...|...|...
492|200|204|205|2|2|3
493|200|205|201|2|2|2
494|200|205|202|2|3|2
495|200|205|203|1|1|1
496|200|205|204|2|1|2

Where 'gid' is group ID, 'aid' is actor ID, 'pid' is partner ID and 'iv', 'm', 'dv' are variables of interest. 

To call the example dataset from mlsrm library, you may use the following code:
```R
data(srm_data)
df <- srm_data
```

If you have your own dataset, you may read the dataset using the following code:

Reading SPSS dataset:
```R
library(foreign)
# For Mac Users
df <- data.frame(read.spss(file.choose()))
# For Windows Users
df <- data.frame(read.spss(choose.files()))
```
Reading CSV dataset:
```R
# For Mac Users
df <- data.frame(read.csv(file.choose()))
# For Windows Users
df <- data.frame(read.csv(choose.files()))
```
## Main Analysis
### (1) Null Model
```R
m0 <- mlsrm(m ~ 1, "gid", "aid", "pid", df)
dv0 <- mlsrm(dv ~ 1, "gid", "aid", "pid", df)
```

### (2) Adding Predictors
```R
m1 <- mlsrm(m ~ iv, "gid", "aid", "pid", df)
dv1 <- mlsrm(dv ~ iv, "gid", "aid", "pid", df)
dv2 <- mlsrm(dv ~ iv + m, "gid", "aid", "pid", df)
```
