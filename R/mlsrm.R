mcmc.onestage <- function(fittedModel, rep=20000){
  suppressWarnings(require(MASS))
  suppressWarnings(require(nlme))
  suppressWarnings(require(lme4))
  stage1_model = fittedModel
  s1.pest = fixed.effects(stage1_model)
  s1.acov = vcov(stage1_model)
  s1.names = names(fixef(stage1_model))
  mcmc <- mvrnorm(rep,s1.pest,s1.acov,empirical=FALSE)
  colnames(mcmc)=s1.names
  return(mcmc)
}

mcmc.twostage <- function(stage1_model,stage2_model, rep=20000){
  suppressWarnings(require(MASS))
  suppressWarnings(require(nlme))
  suppressWarnings(require(lme4))
  s1.pest = fixed.effects(stage1_model)
  s2.pest = fixed.effects(stage2_model)
  s1.acov = vcov(stage1_model)
  s2.acov = vcov(stage2_model)
  s1.names = paste0(names(fixef(stage1_model)), ".s1")
  s2.names = paste0(names(fixef(stage2_model)), ".s2")
  temp1 = matrix(0, length(names(fixef(stage1_model))), length(names(fixef(stage2_model))))
  temp2 = matrix(0, length(names(fixef(stage2_model))), length(names(fixef(stage1_model))))
  s1.acov_m = cbind(s1.acov,temp1)
  s2.acov_m = cbind(temp2,s2.acov)
  acov = rbind(s1.acov_m, s2.acov_m)
  pest = c(s1.pest, s2.pest)
  mcmc <- mvrnorm(rep,pest,acov,empirical=FALSE)
  colnames(mcmc)=c(s1.names, s2.names)
  return(mcmc)
}
mcmc.ci <- function(mcmc, conf=95, title = 'Distribution of Effect'){
  low=(1-conf/100)/2
  upp=((1-conf/100)/2)+(conf/100)
  LL=quantile(mcmc,low)
  UL=quantile(mcmc,upp)
  LL4=format(LL,digits=4)
  UL4=format(UL,digits=4)
  cat("Estimate:\t",mean(mcmc),"\n")
  cat("Std.Error:\t",sd(mcmc),"\n")
  cat("t-value:\t",mean(mcmc)/sd(mcmc),"\n")
  cat("p-value:\t",2*pnorm(-abs(mean(mcmc)/sd(mcmc))),"\n")
  cat(conf,"% CI:\t", "[",LL4, ",",UL4,"]\n")
  hist(mcmc,breaks='FD',col='skyblue',xlab=paste(conf,'% Confidence Interval ','LL',LL4,'  UL',UL4), main=title)
}

srm.modelcompare <- function(model.1, model.2, null = NULL, df = NULL, digits = 3){
  line.1 = c("","Model.1", "Model.2")
  AIC = c(AIC(model.1),AIC(model.2))
  BIC = c(BIC(model.1),BIC(model.2))
  LogLik = c(logLik(model.1), logLik(model.2))
  Deviance = c(-2*logLik(model.1), -2*logLik(model.2))
  DeltaDeviance = c(NA, -2*logLik(model.1)+2*logLik(model.2))

  cat("Note:\tDeviance = -2*logLik.\n")
  result = do.call("rbind",
                   list(AIC=AIC,
                        BIC=BIC,
                        LogLik = LogLik,
                        Deviance = Deviance,
                        Delta.Deviance = DeltaDeviance
                        ))
  if(!is.null(df)){
    Sig = c(NA, 1-pchisq(-2*logLik(model.1)+2*logLik(model.2),df))
    result = rbind(result, Sig=Sig)
    cat("Note:\tSig is obtained by chi-square difference test.\n")
  }

  if(!is.null(null)){
    Pseudo.R2 = c(srm.pseudor2(model.1,null), srm.pseudor2(model.2,null))
    Delta.Pseudo.R2 = c(NA, srm.pseudor2(model.2,null)-srm.pseudor2(model.1,null))
    Effect.Size = c(NA, (srm.pseudor2(model.2,null)-srm.pseudor2(model.1,null))/(1-srm.pseudor2(model.2,null)))

    result = do.call("rbind",
                     list(result,
                          Pseudo.R2=Pseudo.R2,
                          Delta.Pseudo.R2 = Delta.Pseudo.R2,
                          Effect.Size = Effect.Size
                     ))
    cat("Note:\tWe follow Kreft & Leeuw (1998) and Singer (1998) formulae\n\tto calculate pseudo R2.\n")
    cat("Note:\tWe use f2 as the effect size measure.\n")
    cat("Note:\tCohen suggests f2 values of 0.02, 0.15, and 0.35 \n\trepresent small, medium, and large effect sizes.\n\n")
  }
  colnames(result) = c("Model.1", "Model.2")
  print(round(result, digits = digits))
}

mlsrm <- function(formula, gid, aid, pid, data) {
  suppressWarnings(require(nlme))

  data_srm	= suppressWarnings(srm.create.dummies(gid, aid, pid, data, merge.original = TRUE))
  mgs = suppressWarnings(srm.maxgrpsize(gid, aid, pid, data))
  temp1 = c(paste("a",1:mgs, sep=""),paste("p",1:mgs, sep=""))
  temp2 = paste(c("~-1",temp1), collapse = ' + ')
  random_effects = as.formula(temp2)

  names(data_srm)[names(data_srm) == gid] ="SID"
  ctrl <- lmeControl(opt='optim')
  srm =
  lme(formula,
      random = list(
        SID = pdBlocked(list(
          pdIdent(~1),
          pdSRM(random_effects)))),
      correlation=corCompSymm(form=~1 | SID/pdSRM_dyad_id),control=ctrl,
      data=data_srm, na.action=na.omit)
  cat("Done! Congratulations!\n")
  sum=summary(srm)
  var=srm.pct(srm)
  print(sum)
  print(var)
  return(srm)
}

srm.pseudor2 <- function(predict, null){
  predict.variances <- as.numeric(VarCorr(predict)[,1])
  predict.num.mem <- (length(predict.variances)-2)/2
  predict.grp.var <- predict.variances[1]
  predict.act.var <- predict.variances[2]
  predict.part.var <- predict.variances[predict.num.mem+2]
  predict.dyd.var <- predict.variances[length(predict.variances)]

  null.variances <- as.numeric(VarCorr(null)[,1])
  null.num.mem <- (length(null.variances)-2)/2
  null.grp.var <- null.variances[1]
  null.act.var <- null.variances[2]
  null.part.var <- null.variances[null.num.mem+2]
  null.dyd.var <- null.variances[length(null.variances)]

  null.total.var <- null.grp.var + null.act.var + null.part.var + null.dyd.var
  null.act.pct <- null.act.var/null.total.var
  null.part.pct <- null.part.var/null.total.var
  null.grp.pct <- null.grp.var/null.total.var
  null.dyd.pct <- null.dyd.var/null.total.var

  r2.grp <- (null.grp.var - predict.grp.var)/null.grp.var
  r2.act <- (null.act.var - predict.act.var)/null.act.var
  r2.part <- (null.part.var - predict.part.var)/null.part.var
  r2.dyd <- (null.dyd.var - predict.dyd.var)/null.dyd.var

  pseudo.r2 <- r2.grp*null.grp.pct + r2.act*null.act.pct + r2.part*null.part.pct + r2.dyd*null.dyd.pct
  #cat("Pseudo R-squared:\t",pseudo.r2,"\n")
  return(pseudo.r2)
}

srm.maxgrpsize <- function(group.id, act.id, part.id, d, include.self=FALSE, merge.original=FALSE) {

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
  d.dt=aggregate(unique_indiv_id~orig_group_id,all.indivs,length)
  colnames(d.dt)=c("orig_group_id","group_size")
  agg <- data.frame(d.dt)
  max_group_size <- max(agg[,2])
  return(max_group_size)
}

srm.r2 <- function(predict, null){
  predict.variances <- as.numeric(VarCorr(predict)[,1])
  predict.num.mem <- (length(predict.variances)-2)/2
  predict.grp.var <- predict.variances[1]
  predict.act.var <- predict.variances[2]
  predict.part.var <- predict.variances[predict.num.mem+2]
  predict.dyd.var <- predict.variances[length(predict.variances)]

  null.variances <- as.numeric(VarCorr(null)[,1])
  null.num.mem <- (length(null.variances)-2)/2
  null.grp.var <- null.variances[1]
  null.act.var <- null.variances[2]
  null.part.var <- null.variances[null.num.mem+2]
  null.dyd.var <- null.variances[length(null.variances)]

  null.total.var <- null.grp.var + null.act.var + null.part.var + null.dyd.var
  null.act.pct <- null.act.var/null.total.var
  null.part.pct <- null.part.var/null.total.var
  null.grp.pct <- null.grp.var/null.total.var
  null.dyd.pct <- null.dyd.var/null.total.var

  r2.grp <- (null.grp.var - predict.grp.var)/null.grp.var
  r2.act <- (null.act.var - predict.act.var)/null.act.var
  r2.part <- (null.part.var - predict.part.var)/null.part.var
  r2.dyd <- (null.dyd.var - predict.dyd.var)/null.dyd.var

  total.r2 <- r2.grp*null.grp.pct + r2.act*null.act.pct + r2.part*null.part.pct + r2.dyd*null.dyd.pct

  return(total.r2)
}


srm.pseudo.r2 <- function(null.model, predict.model) {

  # Get the variances for null using VarCorr
  variances.null <- as.numeric(VarCorr(null.model)[,1])
  num.mem <- (length(variances.null)-2)/2
  grp.var.null <- variances.null[1]
  act.var.null <- variances.null[2]
  part.var.null <- variances.null[num.mem+2]
  dyd.var.null <- variances.null[length(variances.null)]
  total.var.null <- grp.var.null + act.var.null + part.var.null + dyd.var.null
  null.vals <- c(grp.var.null, act.var.null, part.var.null, dyd.var.null, total.var.null)

  # Get the variances for predict using VarCorr
  variances.predict <- as.numeric(VarCorr(predict.model)[,1])
  num.mem <- (length(variances.predict)-2)/2
  grp.var.predict <- variances.predict[1]
  act.var.predict <- variances.predict[2]
  part.var.predict <- variances.predict[num.mem+2]
  dyd.var.predict <- variances.predict[length(variances.predict)]
  total.var.predict <- grp.var.predict + act.var.predict + part.var.predict + dyd.var.predict
  predict.vals <- c(grp.var.predict, act.var.predict, part.var.predict, dyd.var.predict, total.var.predict)

  # put it together
  tab <- data.frame(null.vals, predict.vals)
  colnames(tab) <- c("null", "predict")

  # Calculate pseudo R sq
  tab$pseudoR2 <- (tab$null-tab$predict)/tab$null

  # Return this
  return(tab)
}