###              A class for the Social Relations Model 
###
### Copyright 2016  Andrew P Knight (knightap@wustl.edu)
### http://apknight.org
### Any errors? Please let me know. With many eyes, all bugs are shallow. 
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

################################################
# This code provides a new method to use with nlme.
# The method enables the estimation of the Kenny's (1994)
# social relations model using the approach described in 
# Snijders & Kenny (1999). This dummy variable approach
# requires placing constraints on the variance-covariance
# matrix.  

# This code was created based on the existing methods
# provided by Pinheiro & Bates in nlme for pdSymm and pdCompSymm

# To run the SRM, you could use the following type of call
# o <- lme(dv ~ 1, random = list(group_id = pdBlocked(list(pdIdent(~1), pdSRM(~a1 + a2 + a3  + a4 + p1 + p2 + p3 + p4-1)))),correlation=corCompSymm(form=~1 | group_id/dyad_id), data=d, na.action=na.omit)
# Where group_id is the grouping identifier; a1-a4 and p1-p4 are the dummies; dyad_id is a dyad marker. Note, following Snijders & Kenny, a1...an and p1...pn, where n = the maximum group size
################################################
require(nlme)

pdSRM <- function (value = numeric(0), form = NULL, nam = NULL, data = sys.frame(sys.parent())) 
{
	object <- numeric(0)
    class(object) <- c("pdSRM", "pdMat")
    pdConstruct(object, value, form, nam, data)
}
environment(pdSRM) <- asNamespace('nlme')

pdConstruct.pdSRM <- function (object, value = numeric(0), form = formula(object), 
    nam = Names(object), data = sys.frame(sys.parent()), ...) 
{
    val <- NextMethod()
    if (length(val) == 0) {
        if ((nc <- length(Names(val))) > 0) {
            attr(val, "ncol") <- nc
        }
        return(val)
    }
    if (is.matrix(val)) {    
    
    	# Read in a Cholesky factorization of the variance-covariance matrix
    	# Then, transform it to the original variance-covariance matrix
        mat.cov <- crossprod(val)
        
        # Build the original correlation matrix based on the variance-covariance matrix
        aux <- 1/sqrt(diag(value))
        mat.cor <- aux * t(mat.cov * aux)
        nc <- dim(mat.cov)[2]
        
        # Extract the variances from the original matrix        
		variances <- diag(mat.cov)		

		# calculate the actor and partner intercepts as the mean of the constituent parts
		a.var <- mean(variances[1:(nc/2)])
		a.sd <- sqrt(a.var)
		p.var <- mean(variances[(nc/2+1):nc])		
		p.sd <- sqrt(p.var)		
				
		# calculate the actor-partner covariance and correlation
		ap.cor <- mean(mat.cor[cbind((1:(nc/2)),(nc/2+1):nc)])        
		ap.cov <- mean(mat.cov[cbind((1:(nc/2)),(nc/2+1):nc)])        

		# Create a new correlation matrix using these parameters
		new.mat.cor <- diag(8)
		new.mat.cor[cbind((1:(nc/2)),(nc/2+1):nc)] <- rep(ap.cor,(nc/2))  	    					
		new.mat.cor[cbind((nc/2+1):nc,(1:(nc/2)))] <- rep(ap.cor,(nc/2))	
		
		# Create the vector of parameters to send to other functions		
		parms <- c(a.sd, p.sd, ap.cor)
		attributes(parms) <- attributes(val)[names(attributes(val)) != "dim"]
        attr(parms, "ncol") <- nc
        class(parms) <- c("pdSRM", "pdMat")        
        return(parms)
    }
}
environment(pdConstruct.pdSRM) <- asNamespace('nlme')

pdMatrix.pdSRM <- function (object, factor = FALSE) 
{
    if (!isInitialized(object)) {
        stop("cannot extract the matrix from an uninitialized \"pdSRM\" object")
    }
    if (is.null(Ncol <- attr(object, "ncol"))) {
        stop("cannot extract the matrix with uninitialized dimensions")
    }
	parms <- as.vector(object)

	# Recreate all the components
	a.sd <- parms[1]
	a.var <- a.sd^2
	p.sd <- parms[2]
	p.var <- p.sd^2
	ap.cor <- parms[3]
	ap.cov <- ap.cor*a.sd*p.sd

	# Create the variance/covariance matrix
	mat.cov <- diag(c(rep(a.var, (Ncol/2)), rep(p.var, (Ncol/2))))
	mat.cov[cbind((1:(Ncol/2)),(Ncol/2+1):Ncol)] <- rep(ap.cov,(Ncol/2))  	    					
	mat.cov[cbind((Ncol/2+1):Ncol,(1:(Ncol/2)))] <- rep(ap.cov,(Ncol/2))		
	if(factor) {
		value <- chol(mat.cov)	
	} else {
		value <- mat.cov
	}
	dimnames(value) <- attr(object, "Dimnames")
	value
}
environment(pdMatrix.pdSRM) <- asNamespace('nlme')


coef.pdSRM <- function (object, unconstrained = TRUE, ...) 
{
  if (unconstrained || !isInitialized(object)) NextMethod()
  else {
    if (is.null(Ncol <- attr(object, "ncol"))) {
      stop("cannot obtain constrained coefficients with uninitialized dimensions")
    }
    val <- as.vector(object)
    val <- c(val[1], val[2], val[3])
    names(val) <- c("std. dev-a","std. dev-p", "corr.")
    val
  }
}
environment(coef.pdSRM) <- asNamespace('nlme')


corMatrix.pdSRM <- function (object, ...) 
{
    if (!isInitialized(object)) {
        stop("cannot extract the matrix from an uninitialized \"pdSRM\" object")
    }
    if (is.null(Ncol <- attr(object, "ncol"))) {
        stop("cannot extract the matrix with uninitialized dimensions")
    }
    obj <- as.vector(object)
    aux <- c(obj[1], obj[2], obj[3])
    
    # This builds the correlation matrix
    value <- diag(Ncol)
	value[cbind((1:(Ncol/2)),(Ncol/2+1):Ncol)] <- rep(aux[3],(Ncol/2))
	value[cbind((Ncol/2+1):Ncol,(1:(Ncol/2)))] <- rep(aux[3],(Ncol/2))
	
	# This builds the vector of SD values    
	attr(value, "stdDev") <- c(rep(aux[1], (Ncol/2)), rep(aux[2], (Ncol/2)))
	attr(value, "corr") <- aux[3]
    if (length(nm <- Names(object)) == 0) {
        nm <- paste("V", 1:Ncol, sep = "")
        dimnames(value) <- list(nm, nm)
    }
    names(attr(value, "stdDev")) <- nm
    value
}
environment(corMatrix.pdSRM) <- asNamespace('nlme')

summary.pdSRM <- function (object, structName = "Social Relations Model", ...) 
{   
    if (isInitialized(object)) {
    	# Build the correlation matrix
        value <- corMatrix(object)
        attr(value, "structName") <- structName
        attr(value, "noCorrelation") <- FALSE
        attr(value, "formula") <- formula(object)
        class(value) <- "summary.pdMat"
        value
    }
    else {
        object
    }    
}
environment(summary.pdSRM) <- asNamespace('nlme')


################################################
# Here is a function for processing the output
# of nlme run using the above method. 
################################################
srm.pct <- function(object) {

	# Get the variances using VarCorr
	variances <- as.numeric(VarCorr(object)[,1])
	num.mem <- (length(variances)-2)/2
	grp.var <- variances[1]
	act.var <- variances[2]
	part.var <- variances[num.mem+2]
	dyd.var <- variances[length(variances)]
	
	# Get the correlations directly from the summary object
	# Haven't written the method for SRM to pull these directly using VarCorr
	# Hack!!
	o.sum <- summary(object)
	
#	ap.cor <- attr(o.sum$apVar, "Pars")[4]
	o <- as.matrix(o.sum$modelStruct$reStruct[[1]])
	ap.cor <- o[(num.mem+2), 2]/sqrt(o[2,2]*o[(num.mem+2),(num.mem+2)])
	ap.cov <- ap.cor*sqrt(act.var*part.var)
	dyd.cor <- coef(object$modelStruct$corStruct,unconstrained=FALSE)

#	dyd.cor <- as.vector(o.sum$modelStruct$corStruct[[1]])/2
	dyd.cov <- dyd.cor*dyd.var
	variance.parms <- as.numeric(c(grp.var, act.var, part.var, dyd.var, ap.cov, dyd.cov))
	names(variance.parms) <- c("Group", "Actor", "Partner", "Dyad", "Generalized Reciprocity", "Dyadic Reciprocity")	
	
	# Compute the percentages and return the correlations
	total.var <- grp.var + act.var + part.var + dyd.var
	act.pct <- 100*act.var/total.var
	part.pct <- 100*part.var/total.var
	grp.pct <- 100*grp.var/total.var
	dyd.pct <- 100*dyd.var/total.var	
	variance.pcts <- c(grp.pct, act.pct, part.pct, dyd.pct, ap.cor, dyd.cor)
	names(variance.pcts) <- c("Group", "Actor", "Partner", "Dyad", "Generalized Reciprocity", "Dyadic Reciprocity")
	output <- as.data.frame(list(variances.and.covariances=variance.parms, percents.and.correlations=variance.pcts))
	return(output)
}
