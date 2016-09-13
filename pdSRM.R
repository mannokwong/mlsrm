###              An nlme class for the Social Relations Model 
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

# Snijders, T. A. B., & Kenny, D. A. 1999. The social relations model for family data: A multilevel approach. Personal Relationships, 6: 471–486.

# This code was created based on the existing methods
# provided by Pinheiro & Bates in nlme for pdSymm and pdCompSymm

# To run the SRM, you could use the following type of call:
# o <- lme(dv ~ 1, random = list(group_id = pdBlocked(list(pdIdent(~1), pdSRM(~a1 + a2 + a3  + a4 + p1 + p2 + p3 + p4-1)))),correlation=corCompSymm(form=~1 | group_id/dyad_id), data=d, na.action=na.omit)
# Where group_id is the grouping identifier; a1-a4 and p1-p4 are the dummies; dyad_id is a dyad marker. Note, following Snijders & Kenny, a1...an and p1...pn, where n = the maximum group size

# Then, you could use the helper function to extract variance parameters and percentages: 
# var.out <- srm.pct(o)
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
        
        # Check to see if this is positive-definite
        
        # Build the original correlation matrix based on the variance-covariance matrix
        aux <- 1/sqrt(diag(mat.cov))
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
		new.mat.cor <- diag(nc)
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

	# Create a correlation matrix
	aux <- 1/sqrt(diag(mat.cov))
	mat.cor <- aux * t(mat.cov * aux)	

	if(factor) {	
		# Test for positive definite here
		cholStatus <- try(u <- chol(mat.cov), silent = TRUE)
		cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE)						
		if(cholError) {
			cat("matrix is not positive definite: executing work around...you should really check your results my friend!\n")
			value <- upper.tri(mat.cov, diag=TRUE)
		} else {
			value <- chol(mat.cov)		
		}
	
		ld <- determinant(mat.cov, logarithm=TRUE)[1]	
		attr(value, "logDet") <- ld$modulus
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
	output <- round(as.data.frame(list(variances.and.covariances=variance.parms, percents.and.correlations=variance.pcts)), 3)
	return(output)
}

################################
# This function creates dummy variables and a unique dyad identifier. The user inputs a 
# dyadic dataset containing at least three id variables: one for group, one for actor, and one for partner. The function will return a dataset containing these original identifiers, plus new ones (using the prefix "pdSRM" to ensure no conflicts with original variable names). The function will also return a set of a... and p... dummies ranging from 1 to n where n is the size of the largest group. 

# Note that this function will return a full dyadic dataset. That is, the function will take note of all individuals who show up either on the actor or the partner side and build a dyadic dataset with a "full" set of dyadic observations for each group. That is, it creates the full round robin, irrespective of whether the original dataset contains these. 

# The arguments for the function are as follows: 

# group.id 	= supply a string with the name of your group identifier variable
# act.id 	= supply a string with the name of your actor identifier variable
# part.id	= supply a string with the name of your partner identifier variable
# include.self 	= logical for whether you want the dataset to contain self ratings. By default the function excludes these. 
# merge.original = logical for whether you want to return the full original dataset, with the new identifiers and dummies added to the dataset or whether you just want a data.frame containing the old and new identifier variables. Default is to just return the identifiers. 

# Sample call is: 

# new.d <- srm.create.dummies(group.id = "team_id", act.id = "act_id", part.id = "part_id", d = old.d)

################################


srm.create.dummies <- function(group.id, act.id, part.id, d, include.self=FALSE, merge.original=FALSE) {
	require(data.table)
	
	# sort the dataset by group.id, act.id, part.id
	d <- d[with(d,order(d[,group.id], d[,act.id], d[, part.id])), ]

	# get just the identifiers
	d.sub <- d[,c(group.id, act.id, part.id)]

	# get the unique groups
	grps <- unique(d.sub[,c(group.id)])
	
	# create a unique actor id and a unique partner id (in case the input data set has these as nested values (rather than unique)
	d.sub$act_indiv_id <- paste(d.sub[,c(group.id)], d.sub[,c(act.id)], sep="_-")
	
	d.sub$part_indiv_id <- paste(d.sub[,c(group.id)], d.sub[,c(part.id)], sep="_-")	
	
	# get the unique individuals as anyone who shows up on the actor or partner side
	acts <- unique(d.sub$act_indiv_id)
	parts <- unique(d.sub$part_indiv_id)	
	all.indivs <- data.frame(unique(c(acts,parts)), stringsAsFactors=F)
	colnames(all.indivs) <- "string_indiv_id"

	# create a unique indiv_id number for everyone
	all.indivs$unique_indiv_id <- 1:length(all.indivs$string_indiv_id)
	
	# Split back apart into the group and indiv identifiers
	all.indivs$orig_group_id <- as.numeric(t(matrix(unlist(strsplit(all.indivs$string_indiv_id, split="_-")), nrow=2))[,1])

	all.indivs$orig_indiv_id <- as.numeric(t(matrix(unlist(strsplit(all.indivs$string_indiv_id, split="_-")), nrow=2))[,2])	
	
	# Using all.indivs and the orig_group_id and the unique_indiv_id create the dyad dataset, with dummies
	
	# get the maximum group size
	d.dt <- data.table(all.indivs)
	agg <- data.frame(d.dt[,list(group_size = length(unique_indiv_id)), by=list(orig_group_id)])	
	max_group_size <- max(agg[,2])
	
	# Create the dyad dataset with an act_number and part_number variable, build it from the ground up
	grps <- unique(all.indivs$orig_group_id)
	grp <- grps[1]
	count <- 1
	for(grp in grps) {
	
		# Get the people in this group
		members <- unique(all.indivs[all.indivs$orig_group_id == grp, c("unique_indiv_id")])
		
		act_num <- 1
		for(act in members) {
		
			part_num <- 1
			for(part in members) {
			
				res.line <- c(grp, act, act_num, part, part_num)
				
				if(count == 1) {
				
					res <- res.line
				
				} else {
				
					res <- rbind(res, res.line)
				
				}
				part_num <- part_num + 1
				count <- count + 1
			
			}
			act_num <- act_num + 1
		
		}		
	
	}
	res <- data.frame(res)
	colnames(res) <- c("orig_group_id", "unique_act_id", "act_num", "unique_part_id", "part_num")

	# create the dummies for each group
	res[,c(paste("a",1:max_group_size, sep=""), paste("p", 1:max_group_size, sep=""))] <- NA
	for(i in 1:max_group_size) {
		v <- paste("a", i, sep="")
		res[,v] <- ifelse(res$act_num == i, 1, 0)

		v <- paste("p", i, sep="")
		res[,v] <- ifelse(res$part_num == i, 1, 0)		
	}
	
	# If self ratings should be excluded...
	if(!include.self) {
	
		res <- res[res$unique_act_id != res$unique_part_id, ]
		
	}
	
	# Add a unique dyad_id
	res$dyad_id <- NA
	count <- 1
	for(grp in grps) {

		# Get the actors and partners in this team
		actors <- unique(res[res$orig_group_id == grp, c("unique_act_id")])
		partners <- unique(res[res$orig_group_id == grp, c("unique_part_id")])	
		indivs <- sort(unique(c(actors, partners)))
		for(a in 1:length(indivs)) {
			for(p in 1:length(indivs)) {
				if(a > p) {
					res$dyad_id <- ifelse((res$unique_act_id == indivs[a] & res$unique_part_id == indivs[p]) | (res$unique_act_id == indivs[p] & res$unique_part_id == indivs[a]), count, res$dyad_id)
					count <- count + 1
				}
			}	
		}
	}	
	
	# Add the original identifiers to this data.frame
	res1 <- merge(res, all.indivs[,c("orig_group_id", "unique_indiv_id", "orig_indiv_id")], by.x=c("orig_group_id", "unique_act_id"), by.y=c("orig_group_id", "unique_indiv_id"), all.x=T)	
	colnames(res1)[length(res1)] <- c("orig_act_id")

	res2 <- merge(res1, all.indivs[,c("orig_group_id", "unique_indiv_id", "orig_indiv_id")], by.x=c("orig_group_id", "unique_part_id"), by.y=c("orig_group_id", "unique_indiv_id"), all.x=T)	
	colnames(res2)[length(res2)] <- c("orig_part_id")	
	
	colnames(res2) <- c(group.id, "pdSRM_part_id", "pdSRM_act_id", "pdSRM_act_num", "pdSRM_part_num",paste("a",1:max_group_size, sep=""), paste("p", 1:max_group_size, sep=""), "pdSRM_dyad_id", act.id, part.id)
	
	# Order the rows columns in a nicer way

	res3 <- res2[with(res2,order(res2[,group.id], res2[,"pdSRM_act_id"], res2[, "pdSRM_part_id"])),c(group.id, act.id, part.id, "pdSRM_act_id", "pdSRM_part_id", "pdSRM_dyad_id", "pdSRM_act_num", "pdSRM_part_num", paste("a",1:max_group_size, sep=""), paste("p", 1:max_group_size, sep=""))]
	
	# Merge with the original dataset
	if(merge.original) {
		res4 <- merge(res3, d, by=c(group.id, act.id, part.id), all.x=T)
		return(res4)
	} else {
		return(res3)
	}
}




