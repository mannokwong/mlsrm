####################################
# pdSRM is the function that gets called by lme. 
# pdSRM is what estimates the random effects
# Begin just by sending an uninitialized object of class pdSRM. 
# Also sends the data from lme
####################################

pdSRM <- function (value = numeric(0), form = NULL, nam = NULL, data = sys.frame(sys.parent())) 
{
	object <- numeric(0)
    class(object) <- c("pdSRM", "pdMat")
    pdConstruct(object, value, form, nam, data)
}
environment(pdSRM) <- asNamespace('nlme')


####################################
# pdConstruct is the function that controls the looping and optimization 
####################################

pdConstruct.pdSRM <- function (object, value = numeric(0), form = formula(object), 
    nam = Names(object), data = sys.frame(sys.parent()), ...) 
{
    val <- NextMethod()
	# Not sure yet what the Next Method does. But, this is where it is getting
	# the input from lme
    
    # If the initial length of val is zero (b/c uninitialized), then 
    # return val with the number of columns in the matrix as an attribute
    # So, this if clause only goes through once - at the very start
    if (length(val) == 0) {
        if ((nc <- length(Names(val))) > 0) {
            attr(val, "ncol") <- nc
        }
        return(val)
    }

    # Once we go through pdConstruct a few times, we end up with a matrix of values
    # We get the Cholesky factorization of the variance-covariance matrix (where is this coming from??)
    if (is.matrix(val)) {    
    	# by computing the crossproduct, of the Cholesky factorization, we get the 
    	# original variance-covariance matrix
    	# crossprod(val) == t(val) %*% val; so, val == chol(value)    	
        value <- crossprod(val)
        
        # Get the dimensionality of the matrix
        nc <- dim(value)[2]
        
        # Convert the matrix into a correlation matrix
        # this is a way to do it that multiplies each
        # cell by the inverse of (a) the standard deviation on one side and (b) the standard deviation on the other side
        aux <- 1/sqrt(diag(value))
        aux <- aux * t(value * aux)
        
        # calculate the average correlation
        aux <- mean(aux[row(aux) != col(aux)])

		# Check to see if the initializing matrix is 
		# positive definite. For it to be positive
		# definite, the average correlation has to be 
		# greater than -1/(nc-1)
		if(aux <= -1/(nc-1)) {      	
        	# If it is not positive definite, force it to be positive definite by constraining the average correlation to be equal to 
            aux <- -1/nc
            warning("initializing \"pdSRM\" object is not positive definite")
        }

		# Set the parameters here as:
		# [1] = the log of the average of the diagonal of the vcov matrix / 2 = log(sd) constrained to equality
		# [2] = the log of a modified fisher transformation of the average correlation in the matrix
        value <- c(log(mean(diag(value)))/2, log((aux + 1/(nc - 1))/(1 - aux)))
        
        # Get the names for the attributes from the inputting object
        attributes(value) <- attributes(val)[names(attributes(val)) != "dim"]
        
        # Set the attribute ncol, which is used below to the number of columns
        attr(value, "ncol") <- nc
        
        # Set the class of this object as a pdSRM, pdMat object
        class(value) <- c("pdSRM", "pdMat")
        
        # Return this pdSRM object, which will now hit the "next method" part and go to build the matrix
        return(value)
    }
}
environment(pdConstruct.pdSRM) <- asNamespace('nlme')

####################################
# pdMatrix is the 
####################################

pdMatrix.pdSRM <- function (object, factor = FALSE) 
{
	# Can't push forward if this is not initialized; probably not an issue
	# since I'm only writing this to use with lme. 
    if (!isInitialized(object)) {
        stop("cannot extract the matrix from an uninitialized \"pdCompSymm\" object")
    }
    if (is.null(Ncol <- attr(object, "ncol"))) {
        stop("cannot extract the matrix with uninitialized dimensions")
    }
    
    # Convert object into a simple vector of values
    obj <- as.vector(object)
    
    # Exponentiate the second element, which is the transformed average correlation
    aux <- exp(obj[2])
    
	# Convert the input obj into a vector of values transformed back
	# Element 1 is the variance, Element 2 is the correlation
    aux <- c(exp(2 * obj[1]), (aux - 1/(Ncol - 1))/(aux + 1))
    
    # Not sure when or why the program dives into the factor clause
    # One place factor is set to "TRUE" is in the logDet function
    if (factor) {
    	# The original code calls the pdFactor function, which does some kind of factoring of the matrix. I don't know what kind. It might be matrix log. I just don't know. The return value from that function needs to be converted from a vector into a matrix. 
    	
    	# What is good news is that I get the same solution if I use a Cholesky factoring of the matrix. So, that's what I've done here. 

# 		 Not using the pdFactor function here; instead doing my own
		
    	# first, create a matrix full of the average correlation
        mat <- array(aux[2], c(Ncol, Ncol))
        # Now make the diagonal a one
        mat[row(mat) == col(mat)] <- 1
        mat <- aux[1] * mat
        
        # This is the variance-covariance matrix
		#value.from.chol <- chol(mat)
		print("mat")
    	print(mat)
		value <- chol(mat)
#        value <- array(pdFactor(object), c(Ncol, Ncol))
		print("value")
		print(value)
		

        # This calculates the logDet; not sure how this is being used
        attr(value, "logDet") <- Ncol * obj[1] + ((Ncol - 1) * 
            log(1 - aux[2]) + log(1 + (Ncol - 1) * aux[2]))/2
    }
    else {
    	# If not factoring it, then just create the matrix of values
    	
    	# first, create a matrix full of the average correlation
        value <- array(aux[2], c(Ncol, Ncol))
        # Now make the diagonal a one
        value[row(value) == col(value)] <- 1
        # Now multiply the variance times the matrix, which gives a variance-covariance matrix
        # Because the correlation times sqrt(varA*varB) (or variance) = covariance
        value <- aux[1] * value
    }
    dimnames(value) <- attr(object, "Dimnames")
    value
}
environment(pdMatrix.pdSRM) <- asNamespace('nlme')

####################################
# This function is called somewhere, but I'm not sure where
# I think this gets called at some point by pdConstruct....
####################################

logDet.pdSRM <- function (object, ...) 
{
    attr(pdMatrix(object, factor = TRUE), "logDet")
}
environment(logDet.pdSRM) <- asNamespace('nlme')


####################################
# This function is called by pdMatrix above. 
# It gets sent the parameters in the bundle (i.e., log(sd) and log(transformed(r))
# It returns a vector of transformed values. I'm not sure what is going on in here.
####################################

pdFactor.pdSRM <- function (object) 
{
 print("nothing to do here now")


    Ncol <- attr(object, "ncol")
    Factor <- double(Ncol*Ncol)
	object <- as.double(object)
	aux <- exp(object[1]) # This is the sd term
	aux1 <- exp(object[2]) # this is the correlation term, fisher transformed
	aux1 <- (aux1-1/(Ncol-1))/(aux1+1) # This is the untransformed correlation

	# Don't know why we're doing this here
	aux2 <- aux * sqrt(1-aux1)
	aux1 = aux * sqrt((1.0 + (Ncol - 1.0) * aux1) / (Ncol))	    

	# I think this is assigning the main diagonal values to slots in the vector
	for(i in 0:(Ncol-1)){
		Factor[(i*Ncol+1)] = aux1
	}
	
	# I think this is assigning the off-diagonal values to the slots they need to be in
	for(i in 1:(Ncol-1)) {
		aux = -aux2/sqrt(i*(i+1))
		for(j in 0:(i-1)) {
			Factor[(i+ (j*Ncol))+1] = aux
		}
		Factor[(i*(Ncol+1))+1] <- -aux*i
	}
	
	# Now, return the vector of values back
	return(Factor)  
    
}
environment(pdFactor.pdSRM) <- asNamespace('nlme')


nlme:::pdFactor.pdSymm


####################################
# Not sure exactly where this is getting called. But it is getting called.
# It's not called in summary, but is getting called above.
# It gets sent the parameters in the bundle (i.e., log(sd) and log(transformed(r))
# It returns a vector of sd and correlation. 
####################################
coef.pdSRM <- function (object, unconstrained = TRUE, ...) 
{
    if (unconstrained || !isInitialized(object)) 
        NextMethod()
    else {
        if (is.null(Ncol <- attr(object, "ncol"))) {
            stop("cannot obtain constrained coefficients with uninitialized dimensions")
        }
                
        val <- as.vector(object)
        aux <- exp(val[2])
        val <- c(exp(val[1]), (aux - 1/(Ncol - 1))/(aux + 1))
        names(val) <- c("std. dev", "corr.")
        val
    }
}
environment(coef.pdSRM) <- asNamespace('nlme')

####################################
# This is used in the summary function.
# It gets sent the parameters in the bundle (i.e., log(sd) and log(transformed(r))
# It returns a vector of the standard deviations and a matrix with the correlations
####################################
corMatrix.pdSRM <- function (object, ...) 
{
    if (!isInitialized(object)) {
        stop("cannot extract the matrix from an uninitialized \"pdSRM\" object")
    }
    if (is.null(Ncol <- attr(object, "ncol"))) {
        stop("cannot extract the matrix with uninitialized dimensions")
    }
    obj <- as.vector(object)
    aux <- exp(obj[2])
    aux <- c(exp(2 * obj[1]), (aux - 1/(Ncol - 1))/(aux + 1))
    value <- array(aux[2], c(Ncol, Ncol))
    value[row(value) == col(value)] <- 1
    attr(value, "stdDev") <- rep(exp(obj[1]), Ncol)
    if (length(nm <- Names(object)) == 0) {
        nm <- paste("V", 1:Ncol, sep = "")
        dimnames(value) <- list(nm, nm)
    }
    names(attr(value, "stdDev")) <- nm
    value
}
environment(corMatrix.pdSRM) <- asNamespace('nlme')

####################################
# Not sure what calls this function; possibly in the summary only
# It gets sent the parameters in the bundle (i.e., log(sd) and log(transformed(r))
# It returns a vector of the standard deviations and a matrix with the correlations
####################################

summary.pdSRM <- function (object, structName = "Social Relations Model (CompSymm template)", ...) 
{
    summary.pdMat(object, structName)
}
environment(summary.pdSRM) <- asNamespace('nlme')