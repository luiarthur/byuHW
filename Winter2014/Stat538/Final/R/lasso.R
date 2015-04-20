source("../../Midterm/PCNSL/R/km.perc.R") #se.pcnt
#rm(list=ls())

library(survival)
library(rms)

library(foreach)
library(doMC)
registerDoMC(16)


##################################################
### read in data

# Clean Data:
  # Expression Data:
  exprss <- read.csv("../Data/BladderCancer_expression.csv")
  # Make the first column the column name
  rownames(exprss) <- exprss[,1]     # These are the gene names
  exprss <- exprss[-nrow(exprss),-1] # Remove the gene names and column of NA's

  # Clinic Data:
  clinic <- read.csv("../Data/BladderCancer_clinical.csv")
  clinic <- clinic[,c("overall.survival","survivalMonth")] # remove patient ID's
  colnames(clinic) <- c("cens","time") 
  clinic$cens <- clinic$cens - 1 # 1 = censored, 2 = death
                                 # Want: 0 = cens, 1 = death

  # No Covariates:
    bigD <- cbind(clinic,t(exprss))
    censor.rate <- mean(bigD$cens == 0)
    N <- nrow(bigD)
    P <- ncol(bigD) - 2

# divide data into training/test sets
  set.seed(1)
  #train.ind <- sample(1:N,round(N/2))
  train.ind <- sample(1:N,130)
  train.set <- bigD[ train.ind,]
  test.set  <- bigD[-train.ind,]
  train.exprss <- train.set[,-c(1:2)]
  test.exprss  <-  test.set[,-c(1:2)]

##################################################
### Exploratory Analysis:
  #KM <- survfit(Surv(train.set$time,train.set$cens)~1,type="kaplan-meier")
  #survplot(KM,xlab="Time (Days)")
  #title("Kaplan Meier Curve")

  one.summary <- function(i) {
    train.summary <- summary(coxph(Surv(train.set$time,train.set$cens) ~ train.exprss[,i]))
    p <- train.summary$coef[5]
    z <- train.summary$coef[4]
    matrix(c("pval"=p,"zscore"=z),1,2)
  }

  system.time(
    pz <- foreach(i=1:ncol(train.exprss),.combine=rbind) %dopar% one.summary(i) # Takes 35s
  )  
  colnames(pz) <- c("pval","zscore")
  pz <- data.frame(pz)
  fdr <- p.adjust(pz$pval,method="fdr")
  
  thresh <- .025 # Only change to get a reasonable number of variables

# identify which markers have fdr less than threshold
  index.nonzero <- which(fdr < thresh)
  length(index.nonzero) # Make sure >= 1?
  gene <- as.matrix(train.exprss[,index.nonzero])

  # fit cox model using selected markers
  marginal.res <- coxph(Surv(train.set$time,train.set$cens) ~ gene)
  summary(marginal.res)
  marginal.res$coef
  #sink("latex/raw/marCox.tex")
  #  xtable(summary(marginal.res)$coef)
  #sink()  
  # Likelihood ratio test=63  on 13 df, p=1.49e-08  n= 130.

  # need to create a coefficient vector (length=number of markers,
  #    all zeros except for selected markers
  marginal.fit.beta <- rep(0,ncol(train.exprss))
  for(i in 1:length(index.nonzero)) {
    marginal.fit.beta[index.nonzero[i]] <- marginal.res$coef[[i]]
  }

  pdf("latex/raw/fdrResid.pdf")
    plot(residuals(marginal.res,type="deviance"),pch=20,col="maroon")
  dev.off()

##################################################
### variable selection under penalized Cox model

  library(penalized)

  # need to format data as dataframe
  train.df <- as.data.frame(train.exprss)
  colnames(train.df) <- as.character(c(1:ncol(train.df)))
  attach(train.df)

  # use cross-validation to obtain penalty parameter (lambda)

  # 5 minutes
  system.time(
    pen.optL1  <- optL1(response=Surv(train.set$time,train.set$cens),
                        penalized=train.df,model="cox",fold=10,standardize=T)
  )

  # 10 minute
  system.time(  
    pen.profL1 <- profL1(response=Surv(train.set$time,train.set$cens),
                         penalized=train.df,model="cox",fold=pen.optL1$fold,
                         standardize=T)
  )

  if (sum(pen.profL1$cvl == max(pen.profL1$cvl)) == 1) {
    lambda1.sel <- pen.profL1$lambda[pen.profL1$cvl==max(pen.profL1$cvl)]
    cvl.sel <- max(pen.profL1$cvl)
  } else {
    lambda1.sel <- pen.profL1$lambda[pen.profL1$cvl==max(pen.profL1$cvl)][1]
    cvl.sel <- max(pen.profL1$cvl)[1]
  }
  #pdf("latex/raw/logliklamb.pdf")
  #  plot(pen.profL1$lambda,pen.profL1$cvl,ylab="log Likelihood",xlab="lambda",
  #       main="Log Likelihood Vs. Lambda",col="orange",pch=20,cex=2)
  #  abline(v=lambda1.sel) #29.0494
  #dev.off()


  pen.res <- penalized(response=Surv(train.set$time,train.set$cens),
                       penalized=train.df,model="cox",lambda1=lambda1.sel,
                       standardize=T)
  pen.res

  pen.fit.beta <- coef(pen.res,"all")
  #pen.fit.beta

### selection path
  pen.res.steps <- penalized(response=Surv(train.set$time,train.set$cens),penalized=train.df,model="cox",lambda1=lambda1.sel,steps=50,standardize=T)

#pdf("latex/raw/selectedVar.pdf")  
#  plotpath(pen.res.steps,lwd=3)
#dev.off()

### check model fit
#residuals(pen.res)
#pdf("latex/raw/lassoResid.pdf")
#  plot(residuals(pen.res),cex=.5,col="purple",main="Residuals Plot")
#dev.off()

# compare coefficients with those that would be obtained from cox model
# (same direction, shrunk toward zero)

gene <- as.matrix(train.exprss[,which(pen.fit.beta!=0)])
lasso.mod <- coxph(Surv(train.set$time,train.set$cens)~gene)
pen.fit.beta[pen.fit.beta != 0]

cn <- colnames(train.exprss[,as.numeric(names(pen.fit.beta[pen.fit.beta != 0]))])
#ILMN_1689037 ILMN_1702933
betahat <- pen.fit.beta[pen.fit.beta != 0]
names(betahat) <- cn
betahat <- t(as.matrix(betahat))

#sink("latex/raw/lassBeta.tex")
#  xtable(betahat)
#sink()

###### ASSESS PERFORMANCE OF TWO APPROACHES: RISK GROUPS ##########################

## risk-groups (check ability of two approaches to more accurately identify
##   high- and low-risk groups
##   make two-group comparisions, check via plots and log-rank test)


####### marginal model (marginal.fit.beta)

# risk score for those in training data set
exp(as.matrix(train.exprss)%*%as.vector(marginal.fit.beta))
# risk score for those in test data set
marginal.risk.score.test <- exp(as.matrix(test.exprss)%*%as.vector(marginal.fit.beta))


# create test set risk groups based on risk score (predicted outcome)
# 0: low-risk group
# 1: high-risk group
risk.group <- ifelse(marginal.risk.score.test < median(marginal.risk.score.test),0,1)

# check whether penalized approach identified markers that
# successfully distinguish between high-risk / low-risk groups
survdiff(Surv(test.set$time,test.set$cens)~c(risk.group))

# plot K-M survival functions for high/low risk
  risk <- ifelse(risk.group==0,"Low","High")
  hl.KM <- survfit(Surv(test.set$time,test.set$cens) ~ risk,type="kaplan-meier")
  se.pcnt(.5,hl.KM) #18.5333

  #pdf("latex/raw/fdrKM.pdf")
  #  survplot(hl.KM,xlab="Time (days)")
  #  title("FDR Model: High-Risk vs. Low-Risk")
  #dev.off()


####### penalized model (pen.fit.beta)

# risk score for those in training data set
  fitted.values(pen.res)
  exp(as.matrix(train.exprss)%*%pen.fit.beta)

# risk score for those in test data set
  pen.risk.score.test <- exp(as.matrix(test.exprss)%*%pen.fit.beta)

# create test set risk groups based on risk score (predicted outcome)
# 0: low-risk group
# 1: high-risk group
  risk.group <- ifelse(pen.risk.score.test < median(pen.risk.score.test),0,1)

# check whether penalized approach identified markers that
# successfully distinguish between high-risk / low-risk groups
  survdiff(Surv(test.set$time,test.set$cens)~c(risk.group))

# plot K-M survival functions for high/low risk based on 3yr survival
  risk <- ifelse(risk.group==0,"Low","High")
  hl.KM.2 <- survfit(Surv(test.set$time,test.set$cens)~risk,type="kaplan-meier")
  

  #pdf("latex/raw/lassoKM.pdf")
  #  survplot(hl.KM.2,xlab="Time (Days)")
  #  title("Penalized Model: High-Risk vs. Low-Risk")
  #dev.off()

###### ASSESS PERFORMANCE OF TWO APPROACHES: TIME-DEPENDENT ROC ##############

library(survivalROC)

#time.seq <- seq(0,16,by=.5)
time.seq <- seq(0,140,by=.5)

marginal.AUC.res <- NA
for(i in 1:length(time.seq)) {
marginal.AUC.res[i] <- survivalROC(Stime=test.set$time,status=test.set$cens,marker=marginal.risk.score.test,predict.time=time.seq[i],method="KM")$AUC
}
marginal.AUC.res
#write.table(marginal.AUC.res,"aucData/marginal.txt",quote=F,col.names=F,row.names=F)

pen.AUC.res <- NA
for(i in 1:length(time.seq)) {
pen.AUC.res[i] <- survivalROC(Stime=test.set$time,status=test.set$cens,marker=pen.risk.score.test,predict.time=time.seq[i],method="KM")$AUC
}
pen.AUC.res
#write.table(pen.AUC.res,"aucData/pen.txt",quote=F,col.names=F,row.names=F)

# I changed this. Is this right?
#pdf("latex/raw/fdrAUC.pdf")
#plot(time.seq,marginal.AUC.res,type="l",ylim=c(0,1),main="FDR Time-Dependent ROC",
#     ylab="AUC",xlab="Time",col="blue",lwd=3)
#dev.off()     

#pdf("latex/raw/lassoAUC.pdf")
#  plot(pen.AUC.res,type="l",col="red",lwd=3)
#dev.off()

#legend("topright",legend=c("lasso","marginal"),col=c("red","blue"),lwd=3)


