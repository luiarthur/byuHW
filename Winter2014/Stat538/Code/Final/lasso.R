library(survival)

##################################################
### read in data

dlbcl.data <- as.matrix(read.table("dlbcl_expressiondata.csv",sep=","))
dlbcl.surv <- as.matrix(read.table("dlbcl_survdata.csv",sep=","))
dim(dlbcl.data)
dlbcl.surv

# divide data into training/test sets
set.seed(538)
train.index <- sample(c(1:240),120)
dlbcl.train.data <- dlbcl.data[,train.index]
dlbcl.train.time <- dlbcl.surv[train.index,1]
dlbcl.train.cens <- dlbcl.surv[train.index,2]

dlbcl.test.data <- dlbcl.data[,-train.index]
dlbcl.test.time <- dlbcl.surv[-train.index,1]
dlbcl.test.cens <- dlbcl.surv[-train.index,2]

##################################################
### exploratory analysis

dlbcl.KM <- summary(survfit(Surv(dlbcl.train.time,dlbcl.train.cens)~1,type="kaplan-meier"))

ymax <- 1
ymin <- 0

dlbcl.KM.time <- c(0,dlbcl.KM$time)
dlbcl.KM.surv <- c(1,dlbcl.KM$surv)

plot(dlbcl.KM.time,dlbcl.KM.surv,ylim=c(ymin,ymax),type="s",lty=1,
xlab="Time",
ylab="Survival",
main="K-M Survival Curve",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)

##################################################
### marker-by-marker analysis

dlbcl.data.pvalue <- dlbcl.data.zscore <- NA
for(i in 1:nrow(dlbcl.train.data)) {
  dlbcl.train.summary <- summary(coxph(Surv(dlbcl.train.time,dlbcl.train.cens) ~ dlbcl.train.data[i,]))
  dlbcl.data.pvalue[i] <-dlbcl.train.summary$coef[5]
  dlbcl.data.zscore[i] <- dlbcl.train.summary$coef[4]
}

dlbcl.data.fdr <- p.adjust(dlbcl.data.pvalue,method="fdr")
thresh <- .15
# identify which markers have fdr less than threshold
which(dlbcl.data.fdr < thresh)
index.nonzero <- which(dlbcl.data.fdr < thresh)
length(index.nonzero)
# fit cox model using selected markers
if(length(index.nonzero) > 1) {
  marginal.res <- coxph(Surv(dlbcl.train.time,dlbcl.train.cens)~t(dlbcl.train.data[index.nonzero,]))
}
if(length(index.nonzero) == 1) {
  marginal.res <- coxph(Surv(dlbcl.train.time,dlbcl.train.cens)~dlbcl.train.data[index.nonzero,])
}
summary(marginal.res)
marginal.res$coef
# need to create a coefficient vector (length=number of markers,
#    all zeros except for selected markers
marginal.fit.beta <- rep(0,nrow(dlbcl.train.data))
for(i in 1:length(index.nonzero)) {
  marginal.fit.beta[index.nonzero[i]] <- marginal.res$coef[[i]]
}



##################################################
### variable selection under penalized Cox model

library(penalized)

# need to format data as dataframe
dlbcl.train.df <- as.data.frame(t(dlbcl.train.data))
colnames(dlbcl.train.df) <- as.character(c(1:ncol(dlbcl.train.df)))
attach(dlbcl.train.df)

# use cross-validation to obtain penalty parameter (lambda)
pen.optL1 <- optL1(response=Surv(dlbcl.train.time,dlbcl.train.cens),penalized=dlbcl.train.df,model="cox",fold=10,standardize=T)
pen.profL1 <- profL1(response=Surv(dlbcl.train.time,dlbcl.train.cens),penalized=dlbcl.train.df,model="cox",fold=pen.optL1$fold,standardize=T)
if (sum(pen.profL1$cvl == max(pen.profL1$cvl)) == 1) {
  lambda1.sel <- pen.profL1$lambda[pen.profL1$cvl==max(pen.profL1$cvl)]
  cvl.sel <- max(pen.profL1$cvl)
} else {
  lambda1.sel <- pen.profL1$lambda[pen.profL1$cvl==max(pen.profL1$cvl)][1]
  cvl.sel <- max(pen.profL1$cvl)[1]
}
plot(pen.profL1$lambda,pen.profL1$cvl)
lambda1.sel

pen.res <- penalized(response=Surv(dlbcl.train.time,dlbcl.train.cens),penalized=dlbcl.train.df,model="cox",lambda1=lambda1.sel,standardize=T)
pen.res

pen.fit.beta <- coef(pen.res,"all")
pen.fit.beta

### selection path
pen.res.steps <- penalized(response=Surv(dlbcl.train.time,dlbcl.train.cens),penalized=dlbcl.train.df,model="cox",lambda1=lambda1.sel,steps=50,standardize=T)

#postscript("dlbcl.path.ps",height=10,width=10)
plotpath(pen.res.steps)
#dev.off()

### check model fit
residuals(pen.res)
plot(residuals(pen.res))

# compare coefficients with those that would be obtained from cox model
# (same direction, shrunk toward zero)
if(sum(pen.fit.beta != 0) > 1) {
  coxph(Surv(dlbcl.train.time,dlbcl.train.cens)~t(dlbcl.train.data[which(pen.fit.beta!=0),]))
}
if(sum(pen.fit.beta != 0) == 1) {
  coxph(Surv(dlbcl.train.time,dlbcl.train.cens)~dlbcl.train.data[which(pen.fit.beta!=0),])
}
pen.fit.beta[pen.fit.beta != 0]




###### ASSESS PERFORMANCE OF TWO APPROACHES: RISK GROUPS ##########################

## risk-groups (check ability of two approaches to more accurately identify
##   high- and low-risk groups
##   make two-group comparisions, check via plots and log-rank test)


####### marginal model (marginal.fit.beta)

# risk score for those in training data set
exp(t(marginal.fit.beta)%*%dlbcl.train.data)
# risk score for those in test data set
exp(t(marginal.fit.beta)%*%dlbcl.test.data)
marginal.risk.score.test <- exp(t(marginal.fit.beta)%*%dlbcl.test.data)

# create test set risk groups based on risk score (predicted outcome)
# 0: low-risk group
# 1: high-risk group
risk.group <- ifelse(marginal.risk.score.test < median(marginal.risk.score.test),0,1)

# check whether penalized approach identified markers that
# successfully distinguish between high-risk / low-risk groups
survdiff(Surv(dlbcl.test.time,dlbcl.test.cens)~c(risk.group))

# plot K-M survival functions for high/low risk
g1.time <- dlbcl.test.time[risk.group==0] # low-risk group
g1.cens <- dlbcl.test.cens[risk.group==0]

g2.time <- dlbcl.test.time[risk.group==1] # high-risk group
g2.cens <- dlbcl.test.cens[risk.group==1]

g1.KM <- survfit(Surv(g1.time,g1.cens)~1,type="kaplan-meier")
g2.KM <- survfit(Surv(g2.time,g2.cens)~1,type="kaplan-meier")
summary(g1.KM)
summary(g2.KM)

ymax <- 1
ymin <- 0

g1.KM.time <- c(0,g1.KM$time)
g1.KM.surv <- c(1,g1.KM$surv)
g2.KM.time <- c(0,g2.KM$time)
g2.KM.surv <- c(1,g2.KM$surv)

g1.KM.upper <- c(1,g1.KM$upper)
g1.KM.lower <- c(1,g1.KM$lower)
g2.KM.upper <- c(1,g2.KM$upper)
g2.KM.lower <- c(1,g2.KM$lower)

plot(g1.KM.time,g1.KM.surv,ylim=c(ymin,ymax),type="s",lty=1,
     xlab="Time",
     ylab="Survival",
     main="Marginal Model: High-Risk vs. Low-Risk",
     xlim=c(0,max(dlbcl.test.time)))
#lines(g1.KM.time,g1.KM.lower,type="s",lty=2)
#lines(g1.KM.time,g1.KM.upper,type="s",lty=2)
lines(g2.KM.time,g2.KM.surv,lty=1,type="s",col=2)
#lines(g2.KM.time,g2.KM.lower,type="s",lty=2,col=2)
#lines(g2.KM.time,g2.KM.upper,type="s",lty=2,col=2)



####### penalized model (pen.fit.beta)

# risk score for those in training data set
fitted.values(pen.res)
exp(t(pen.fit.beta)%*%dlbcl.train.data)
# risk score for those in test data set
exp(t(pen.fit.beta)%*%dlbcl.test.data)
pen.risk.score.test <- exp(t(pen.fit.beta)%*%dlbcl.test.data)

# create test set risk groups based on risk score (predicted outcome)
# 0: low-risk group
# 1: high-risk group
risk.group <- ifelse(pen.risk.score.test < median(pen.risk.score.test),0,1)

# check whether penalized approach identified markers that
# successfully distinguish between high-risk / low-risk groups
survdiff(Surv(dlbcl.test.time,dlbcl.test.cens)~c(risk.group))

# plot K-M survival functions for high/low risk based on 3yr survival
g1.time <- dlbcl.test.time[risk.group==0] # low-risk group
g1.cens <- dlbcl.test.cens[risk.group==0]

g2.time <- dlbcl.test.time[risk.group==1] # high-risk group
g2.cens <- dlbcl.test.cens[risk.group==1]

g1.KM <- survfit(Surv(g1.time,g1.cens)~1,type="kaplan-meier")
g2.KM <- survfit(Surv(g2.time,g2.cens)~1,type="kaplan-meier")
summary(g1.KM)
summary(g2.KM)

ymax <- 1
ymin <- 0

g1.KM.time <- c(0,g1.KM$time)
g1.KM.surv <- c(1,g1.KM$surv)
g2.KM.time <- c(0,g2.KM$time)
g2.KM.surv <- c(1,g2.KM$surv)

g1.KM.upper <- c(1,g1.KM$upper)
g1.KM.lower <- c(1,g1.KM$lower)
g2.KM.upper <- c(1,g2.KM$upper)
g2.KM.lower <- c(1,g2.KM$lower)

plot(g1.KM.time,g1.KM.surv,ylim=c(ymin,ymax),type="s",lty=1,
     xlab="Time",
     ylab="Survival",
     main="Penalized Model: High-Risk vs. Low-Risk",
     xlim=c(0,max(dlbcl.test.time)))
#lines(g1.KM.time,g1.KM.lower,type="s",lty=2)
#lines(g1.KM.time,g1.KM.upper,type="s",lty=2)
lines(g2.KM.time,g2.KM.surv,lty=1,type="s",col=2)
#lines(g2.KM.time,g2.KM.lower,type="s",lty=2,col=2)
#lines(g2.KM.time,g2.KM.upper,type="s",lty=2,col=2)




###### ASSESS PERFORMANCE OF TWO APPROACHES: TIME-DEPENDENT ROC ##########################

library(survivalROC)

time.seq <- seq(0,16,by=.5)

marginal.AUC.res <- NA
for(i in 1:length(time.seq)) {
marginal.AUC.res[i] <- survivalROC(Stime=dlbcl.test.time,status=dlbcl.test.cens,marker=marginal.risk.score.test,predict.time=time.seq[i],method="KM")$AUC
}
marginal.AUC.res

pen.AUC.res <- NA
for(i in 1:length(time.seq)) {
pen.AUC.res[i] <- survivalROC(Stime=dlbcl.test.time,status=dlbcl.test.cens,marker=pen.risk.score.test,predict.time=time.seq[i],method="KM")$AUC
}
pen.AUC.res

plot(marginal.AUC.res,type="l",
     ylim=c(0,1),
     main="Time-Dependent ROC",
     ylab="AUC",
     xlab="Time")
lines(pen.AUC.res,type="l",col="red")





