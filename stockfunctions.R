# Setup the environment
apk.r.setup <- function() {
	pckgs <- c("lme4", "psych", "quantreg", "xtable", "reshape", "arm", "MCMCpack", "multilevel", "sqldf", "RMySQL", "lavaan", "Hmisc", "sna", "R2wd")
	install.packages(pckgs, dependencies=TRUE, repos="http://cran.wustl.edu")
}
# Frequently used libraries
require(lme4)
require(psych)
require(quantreg)
require(xtable)
require(reshape)
require(arm)
require(MCMCpack)
require(multilevel)
require(sqldf)
require(metafor)
require(lavaan)
require(sna)
# require(R2wd)

# FUNCTION FOR OUTPUTTING A CORRELATION MATRIX
correlation.matrix <- function(x, diagval=NA){
require(Hmisc)
x <- as.matrix(x)
#j <- as.matrix(j)
R <- rcorr(x)$r
p <- rcorr(x)$P

M <- formatC(round(colMeans(x, na.rm=T),2), digits=2, format="f")
SD <- formatC(round(sd(x, na.rm=T),2),digits=2, format="f")

## define notations for significance levels; spacing is important.
mystars <- ifelse(p < .01, "**", ifelse(p < .05, "* ", ifelse(p < .10, "+ ", " ")))

## trunctuate the matrix that holds the correlations to two decimal
R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]

## build a new matrix that includes the correlations with their apropriate stars
Rnew <- R
Rnew<-matrix(paste(Rnew, mystars, sep=""), ncol=ncol(x))
diag(Rnew) <- diagval
rownames(Rnew) <- colnames(x)
colnames(Rnew) <- paste(colnames(x), "", sep="")

## remove upper triangle
Rnew <- as.matrix(Rnew)
Rnew[upper.tri(Rnew, diag = FALSE)] <- ""

## Add the mean and the SD
Rnew0<-cbind(M, SD, Rnew)

Rnew1 <- as.data.frame(Rnew0)

## remove last column and return the matrix (which is now a data frame)
#Rnew <- cbind(Rnew[1:length(Rnew)-1])
return(Rnew1)
} 

plot.interaction <- function(dv, v1, v2, controls,dfile){
eq<-paste(dv,v1,sep="~")
eq<-paste(eq,v2,sep="*")
if(length(controls) > 0)
{
for(i in 1:length(controls))
{
	eq<-paste(eq,controls[i],sep="+")	
}
}
arg<-model.frame(as.formula(eq), dfile)
out<-glm.nb(arg)

v1.x<-mean(dfile[,c(v1)], na.rm=TRUE)
v1.sd<-sd(dfile[,c(v1)], na.rm=TRUE)
v1.lo<-v1.x-v1.sd
v1.hi<-v1.x+v1.sd

v2.x<-mean(dfile[,c(v2)], na.rm=TRUE)
v2.sd<-sd(dfile[,c(v2)], na.rm=TRUE)
v2.lo<-v2.x-v2.sd
v2.hi<-v2.x+v2.sd

lolo<-out$coefficients[1]+out$coefficients[2]*v1.lo+out$coefficients[3]*v2.lo+out$coefficients[length(out$coefficients)]*v1.lo*v2.lo
lohi<-out$coefficients[1]+out$coefficients[2]*v1.lo+out$coefficients[3]*v2.hi+out$coefficients[length(out$coefficients)]*v1.lo*v2.hi
hilo<-out$coefficients[1]+out$coefficients[2]*v1.hi+out$coefficients[3]*v2.lo+out$coefficients[length(out$coefficients)]*v1.hi*v2.lo
hihi<-out$coefficients[1]+out$coefficients[2]*v1.hi+out$coefficients[3]*v2.hi+out$coefficients[length(out$coefficients)]*v1.hi*v2.hi

loline<-c(lolo, lohi)
hiline<-c(hilo, hihi)

#quartz(type="pdf", file=filename)
plot(hiline, type="o",col="red", ylim=c(1,2.5), axes=FALSE, main="", cex.main=.75, xlab="# of Dispensed", pch=16, ylab="Likelihood of Error Report", lty=1, lwd=2)
axis(1, at=1:2, lab=c("Minus 1 SD", "Plus 1 SD"), tck=0, cex.axis=.6)
axis(2, at=c("2.25", "2.75", "3.25"), lab=c("2.25", "2.75", "3.25"), tck=0, cex.axis=.7)
lines(loline, type="o", col="blue", pch=18, lwd=2, lty=2)
legend("bottomright", legend=c("+1 SD Number of People", "-1 SD Number of People"), cex=.8,col=c("red", "blue"), pch=c(16,18), lty=c(1,2), bty="n")
#dev.off()
return(out)
}

blau<-function(catvar, groupid, dfile){
require(sqldf)
q0<-paste("SELECT ",groupid,", COUNT(",groupid,") as groupsize FROM dfile WHERE ",catvar," IS NOT NULL GROUP BY ",groupid, sep=" ")
out0<-data.frame(sqldf(q0))
new<-merge(dfile, out0, by=c(groupid))
q1<-paste("SELECT ",groupid,",groupsize, COUNT(", catvar,") AS catsize FROM new GROUP BY ",groupid,", ",catvar, sep=" ")
out1<-data.frame(sqldf(q1))
out1$prop2<-(out1$catsize/out1$groupsize)*(out1$catsize/out1$groupsize)
q2<-paste("SELECT ",groupid,", (1-SUM(prop2)) AS ",catvar,"_blau FROM out1 GROUP BY ",groupid,sep="")
d<-data.frame(sqldf(q2))
return<-d[,c(groupid,paste(catvar,"_blau",sep=""))]
}

regtable<-function(models, rnames)
{
	numrows<-length(rnames)	
	rnames<-data.frame(rnames)
	colnames(rnames)<-c("Var")
	df<-c("DF",rep("",length(models)))
    rsq<-c("R2",rep("",length(models)))
    fstat<-c("F",rep("",length(models)))	
	for(i in 1:length(models))
	{	
		s<-mapply(summary, models[i], USE.NAMES = TRUE)
		coef<-s[4,1]$coefficients
		column<-rep("",length(coef[,1]))
		name<-rep("",length(coef[,1]))
		for(z in 1:length(coef[,1]))
		{
			name[z]<-rownames(coef)[z]		
			star <- ifelse(coef[z,4] < .01, "**", ifelse(coef[z,4] < .05, "*", ifelse(coef[z,4] < 0.10, "+","")))
   			val<-formatC(round(coef[z,1],2),digits=2, format="f")
			column[z]<-paste(val,star,sep="")
		}
		coefs<-data.frame(cbind(name, column))
		colnames(coefs)<-c("Var",paste("M",i,sep=""))
		
		df[(i+1)]<-paste(s[10,1]$fstatistic[2], s[10,1]$fstatistic[3], sep=",")
		rsq[(i+1)]<-formatC(round(s[8,1]$r.squared[1],2), digits=2, format="f")
		fval<-formatC(round(s[10,1]$fstatistic[1],2), digits=2, format="f")
		pval<-pf(s[10,1]$fstatistic[1], s[10,1]$fstatistic[2], s[10,1]$fstatistic[3], lower.tail=FALSE)
		star <- ifelse(pval < .001, "***", ifelse(pval < .01, "** ", ifelse(pval < .05, "* ", " ")))
		fstat[(i+1)]<-paste(fval, star, sep="")		
		if(i == 1)
		{
			tab<-merge(coefs, rnames, by=c("Var"), all.y=TRUE, sort=FALSE)
		} else {
			tab<-merge(tab, coefs, by=c("Var"), all.y=TRUE, sort=FALSE)
		}
	}
	
	# Add the summary statistics for each model
	
	data.frame(tab)
	# Now make sure that the order is as desired
	tab<-merge(rnames, tab, by=c("Var"), all.x=TRUE, sort=FALSE)
	
	sumstats<-data.frame(rbind(df, rsq, fstat))
#	sumstats<-sumstats[,c(1,rev(2:length(colnames(tab))))]
	colnames(sumstats)<-colnames(tab)
	tab<-rbind(tab,sumstats)
	return(tab)
}

hlmtable<-function(models, rnames)
{
	numrows<-length(rnames)	
	rnames<-data.frame(rnames)
	colnames(rnames)<-c("Var")
	df<-c("DF",rep("",length(models)))
    rsq<-c("R2",rep("",length(models)))
    fstat<-c("F",rep("",length(models)))	
	for(i in 1:length(models))
	{	
		s<-mapply(summary, models[i], USE.NAMES = TRUE)
		coef<-s[4,1]$coefficients
		column<-rep("",length(coef[,1]))
		name<-rep("",length(coef[,1]))
		for(z in 1:length(coef[,1]))
		{
			name[z]<-rownames(coef)[z]		
			star <- ifelse(coef[z,4] < .01, "**", ifelse(coef[z,4] < .05, "*", ifelse(coef[z,4] < 0.10, "+","")))
   			val<-formatC(round(coef[z,1],2),digits=2, format="f")
			column[z]<-paste(val,star,sep="")
		}
		coefs<-data.frame(cbind(name, column))
		colnames(coefs)<-c("Var",paste("M",i,sep=""))
		
		df[(i+1)]<-paste(s[10,1]$fstatistic[2], s[10,1]$fstatistic[3], sep=",")
		rsq[(i+1)]<-formatC(round(s[8,1]$r.squared[1],2), digits=2, format="f")
		fval<-formatC(round(s[10,1]$fstatistic[1],2), digits=2, format="f")
		pval<-pf(s[10,1]$fstatistic[1], s[10,1]$fstatistic[2], s[10,1]$fstatistic[3], lower.tail=FALSE)
		star <- ifelse(pval < .001, "***", ifelse(pval < .01, "** ", ifelse(pval < .05, "* ", " ")))
		fstat[(i+1)]<-paste(fval, star, sep="")		
		if(i == 1)
		{
			tab<-merge(coefs, rnames, by=c("Var"), all.y=TRUE, sort=FALSE)
		} else {
			tab<-merge(tab, coefs, by=c("Var"), all.y=TRUE, sort=FALSE)
		}
	}
	
	# Add the summary statistics for each model
	
	data.frame(tab)
	# Now make sure that the order is as desired
	tab<-merge(rnames, tab, by=c("Var"), all.x=TRUE, sort=FALSE)
	
	sumstats<-data.frame(rbind(df, rsq, fstat))
#	sumstats<-sumstats[,c(1,rev(2:length(colnames(tab))))]
	colnames(sumstats)<-colnames(tab)
	tab<-rbind(tab,sumstats)
	return(tab)
}

# Output a nice table of results for bayesian models

bayestable <- function(vars) {
	p<-matrix(rep(NA,(6*ncol(vars))), nrow=ncol(vars), ncol=6)
	for(i in 1:ncol(vars)){
		p[i,1] <- mean(vars[,i], na.rm=TRUE)
		p[i,2] <- sd(vars[,i], na.rm=TRUE)
		p[i,4] <- quantile(vars[,i], .025)
		p[i,5] <- quantile(vars[,i], .975)
		if(p[i,1] > 0)
		{
			num <- length(which(vars[,i] < 0))
		}
		else
		{
			num <- length(which(vars[,i] > 0))
		}
		p[i,6] <- num/length(vars[,i])
	}	
	p[,3] <- p[,1] / p[,2]
	rownames(p) <- colnames(vars)
	colnames(p) <- c("X","SD","t",".025",".975","p")
	return(p)
}

# A function for calculating simple slopes given an object output from lme, z is a vector of three points at which you want slopes


# A function for calculating simple slopes given coefficients, standard errors, and the positions in the array of the relevant variables

ss.lin.ols <- function(coefs, ses, x.n, z.n, xz.n,z) {
	b.1 <- coefs[x.n]
	b.2 <- coefs[z.n]
	b.3 <- coefs[xz.n]
	
	z.hi <- mean(z, na.rm=TRUE) + sd(z, na.rm=TRUE)
	z.av <- mean(z, na.rm=TRUE)
	z.lo <- mean(z, na.rm=TRUE) - sd(z, na.rm=TRUE)	
	
	ss.hi <- (b.1+b.3*z.hi)
	ss.av <- (b.1+b.3*z.av)
	ss.lo <- (b.1+b.3*z.lo)
	
	se.hi <- sqrt(ses[x.n,x.n] + 2*z.hi*ses[x.n,xz.n] + (z.hi^2)*ses[xz.n,xz.n])
	se.av <- sqrt(ses[x.n,x.n] + 2*z.av*ses[x.n,xz.n] + (z.av^2)*ses[xz.n,xz.n])
	se.lo <- sqrt(ses[x.n,x.n] + 2*z.lo*ses[x.n,xz.n] + (z.lo^2)*ses[xz.n,xz.n])
	
	tval.hi <- ss.hi/se.hi
	tval.av <- ss.av/se.av 
	tval.lo <- ss.lo/se.lo	
	
	dfs <- length(z) - (length(coefs)-1) - 1
	
	pval.hi <- 2*pt(-abs(tval.hi),df=dfs)
	pval.av <- 2*pt(-abs(tval.av),df=dfs)
	pval.lo <- 2*pt(-abs(tval.lo),df=dfs)	
	
	ss <- rbind(ss.hi, ss.av, ss.lo)
	se <- rbind(se.hi, se.av, se.lo)
	tval <- rbind(tval.hi, tval.av, tval.lo)
	dfs <- rbind(dfs, dfs, dfs)
	pval <- rbind(pval.hi, pval.av, pval.lo)
	
	ret <- cbind(ss, se, tval, dfs, pval)
	row.names(ret) <- c("+1SD Z","Mean Z","-1SD Z")
	colnames(ret) <- c("SS","SE","t","df","p")
	return(ret)
}

ss.lin.rcm <- function(obj, x.n, z.n, xz.n, z) {
	coefs <- obj$coefficients$fixed
	ses <- obj$varFix
	dfs <- obj$fixDF$X
	
	b.1 <- coefs[x.n]
	b.2 <- coefs[z.n]
	b.3 <- coefs[xz.n]
	
	z.hi <- mean(z, na.rm=TRUE) + sd(z, na.rm=TRUE)
	z.av <- mean(z, na.rm=TRUE)
	z.lo <- mean(z, na.rm=TRUE) - sd(z, na.rm=TRUE)		
	
	ss.hi <- (b.1+b.3*z.hi)
	ss.av <- (b.1+b.3*z.av)
	ss.lo <- (b.1+b.3*z.lo)
	
	se.hi <- sqrt(ses[x.n,x.n] + 2*z.hi*ses[x.n,xz.n] + (z.hi^2)*ses[xz.n,xz.n])
	se.av <- sqrt(ses[x.n,x.n] + 2*z.av*ses[x.n,xz.n] + (z.av^2)*ses[xz.n,xz.n])
	se.lo <- sqrt(ses[x.n,x.n] + 2*z.lo*ses[x.n,xz.n] + (z.lo^2)*ses[xz.n,xz.n])
	
	tval.hi <- ss.hi/se.hi
	tval.av <- ss.av/se.av 
	tval.lo <- ss.lo/se.lo	
	
	dfs <- dfs[xz.n]
	
	pval.hi <- 2*pt(-abs(tval.hi),df=dfs)
	pval.av <- 2*pt(-abs(tval.av),df=dfs)
	pval.lo <- 2*pt(-abs(tval.lo),df=dfs)	
	
	ss <- rbind(ss.hi, ss.av, ss.lo)
	se <- rbind(se.hi, se.av, se.lo)
	tval <- rbind(tval.hi, tval.av, tval.lo)
	dfs <- rbind(dfs, dfs, dfs)
	pval <- rbind(pval.hi, pval.av, pval.lo)
	
	ret <- cbind(ss, se, tval, dfs, pval)
	row.names(ret) <- c("+1SD Z","Mean Z","-1SD Z")
	colnames(ret) <- c("SS","SE","t","df","p")
	return(ret)
}

# Function to return typical scale statistics for group-level constructs
scale.scores <- function(vars, dat) {
	count <- 1
	for(i in vars) {
		val <- data.frame(rowMeans(dat[,c(i[[3]])], na.rm=TRUE))
		colnames(val) <- i[[2]]
		if(count == 1) {
			return.dat <- val	
		} else {
			return.dat <- cbind(return.dat, val)
		}
		count <- count+1
	}
	dat <- cbind(dat, return.dat)
	return(dat)
}

scale.stats <- function (vars, groupvar, dat) {
	icc1 <- rep(NA, length(vars))
	icc2 <- rep(NA, length(vars))
	alpha <- rep(NA, length(vars))
	rwgjmed <- rep(NA, length(vars))
	rwgjx <- rep(NA, length(vars))
	rwgjsd <- rep(NA, length(vars))
	rwgjmin <- rep(NA, length(vars))
	rwgjmax <- rep(NA, length(vars))
	name <- rep(NA, length(vars))
	count <- 1
	for(i in vars) {
			name[count] <- i[[1]]
			dv <- dat[, c(i[[2]])]
			gv <- dat[, c(groupvar)]    
			items <- dat[,c(i[[3]])]
			out <- aov(dv ~ as.factor(gv), data=dat)
			icc1[count] <-ICC1(out)
			icc2[count] <-ICC2(out)
			if(length(i[[3]]) <= 1) {
				alpha[count] <- NA
				out<-rwg(items,gv, i[[4]])
				rwgjmed[count]<-median(out$rwg, na.rm=TRUE)			
				rwgjx[count]<-mean(out$rwg, na.rm=TRUE)
				rwgjsd[count]<-sd(out$rwg, na.rm=TRUE)
				rwgjmin[count]<-min(out$rwg, na.rm=TRUE)
				rwgjmax[count]<-max(out$rwg, na.rm=TRUE)        								
			} else {
				out <- alpha(items, na.rm=TRUE, check.keys=FALSE)
				alpha[count] <- out$total['raw_alpha']   		
				out<-rwg.j(items,gv, i[[4]])
				rwgjmed[count]<-median(out$rwg.j, na.rm=TRUE)			
				rwgjx[count]<-mean(out$rwg.j, na.rm=TRUE)
				rwgjsd[count]<-sd(out$rwg.j, na.rm=TRUE)
				rwgjmin[count]<-min(out$rwg.j, na.rm=TRUE)
				rwgjmax[count]<-max(out$rwg.j, na.rm=TRUE)        				
			}			
			count <- count + 1
	}
	v <- cbind(name,  alpha, icc1, icc2, rwgjx, rwgjsd, rwgjmed, rwgjmin, rwgjmax)
	return(v)
}

# This function takes the construct listing and generates two tables - one containing fit indices, the second containing standardized loadings
scale.cfa <- function(vars,dat) {
	count <- 1
	for(construct in vars) {
		itemlist <- paste0(construct$items, collapse=" + ")
		model <- paste('f1 =~ ', itemlist, sep="")
		mod <- cfa(model, data=dat)
		mod.fit <- fitMeasures(mod)
		mod.est <- standardizedSolution(mod)
		fit.ind <- c(round(mod.fit[c(1,2,3,7,8,12,14,16,19,20)],3))
		ind.names <- mod.est[1:length(construct$items),3]
		std.loads <- round(mod.est[1:length(construct$items),4],3)
		var.name <- c(construct$name)
		lds <- cbind(rep(var.name), length(ind.names), ind.names, std.loads)
		if(count == 1) {
			fit.out <- c(var.name, fit.ind)
			load.out <- lds
		} else {
			fit.out <- rbind(fit.out, c(var.name, fit.ind))
			load.out <- rbind(load.out, lds)
		}
		count <- count + 1
	}
	row.names(fit.out) <- 1:length(fit.out[,1])
	fit.out <- data.frame(fit.out)
	load.out <- data.frame(load.out)
	return(list(fit.out, load.out))
}


# Aggregation Function

ak_aggregate <- function(dat, constructs, groupvars) {

	aggvars <- rep(NA, length(constructs))
	for(i in 1:length(constructs)) {
		aggvars[i] <- constructs[[i]]$var	
	}	
	group_list <- vector("list", length(groupvars))
	for(i in 1:length(groupvars)) {
		group_list[[i]] <- dat[,groupvars[i]]
	}
	
	q <- paste("SELECT ",paste(groupvars, collapse=", "), ", COUNT(", groupvars[length(groupvars)],") AS groupsize FROM dat GROUP BY ", paste(groupvars,collapse=","))

	n <- sqldf(q)
	
	agg.x <- aggregate(dat[,aggvars], group_list, mean, na.rm=TRUE)
	colnames(agg.x) <- c(groupvars, paste(aggvars, "_x", sep=""))
	
	agg.sd <- aggregate(dat[,aggvars], group_list, sd, na.rm=TRUE)
	colnames(agg.sd) <- c(groupvars, paste(aggvars, "_sd", sep=""))	
	
	agg.dat <- merge(agg.x, agg.sd, by=groupvars)	
	agg.dat <- merge(agg.dat, n, by=groupvars)
	return(agg.dat)
}
