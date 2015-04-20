library(survival)
library(superpc)
source("../../Midterm/PCNSL/R/km.perc.R") #se.pcnt

##################################################
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

### can't do PCA when p > n: prcomp just ignores any p that are beyond n
### stdev matrix should be a pxp matrix 
###???
#res.prcomp <- prcomp(t(train.exprss),center=T,scale=T,retx=T)
#res.prcomp$sdev
#
#res.prcomp <- prcomp(train.exprss,center=T,scale=T,retx=T)
#res.prcomp$sdev

###############################################
#### superpc ###############
featurenames <- paste("feature",as.character(1:N),sep="")

# create train and test data objects. censoring.status=1 means the event occurred;
# censoring.status=0 means censored

train.pc <- list(x=t(train.exprss),y=train.set$time,censoring.status=train.set$cens,
                 featurenames=featurenames)
test.pc  <- list(x=t(test.exprss),y=test.set$time,censoring.status=test.set$cens,
                 featurenames= featurenames)

# train  the model. This step just computes the scores for each feature
train.obj <- superpc.train(train.pc, type="survival")

# note for regression (non-survival) data, we leave the component "censoring.status"
# out of the data object, and call superpc.train with type="regression".
# otherwise the superpc commands are all the same


# cross-validate the model (i.e., identify Cox score threshold at which 
#                                 markers should be retained)
system.time( 
  cv.obj <- superpc.cv(train.obj, train.pc) # Takes about 1 minute
)

# plot the cross-validation curves. From this plot we see that the 1st 
# principal component is significant and the best threshold is around 1.5
#pdf("latex/raw/superpc.pdf")
#  superpc.plotcv(cv.obj) # 2.5?
#dev.off()

# here we have the luxury of  test data, so we can compute the likelihood ratio statistic
# over the test data and plot them. We see that the threshold of 1.5
# works pretty well (see JASA paper, e.g., section 5)
# This is the same as before but uses Cross Validation
lrtest.obj <- superpc.lrtest.curv(train.obj, train.pc, test.pc)
#pdf("latex/raw/superpcLR.pdf")
#  superpc.plot.lrtest(lrtest.obj)
#dev.off()  
thresh <- lrtest.obj$threshold[which.max(lrtest.obj$lrtest)] # 1.34

# now we derive the predictor of survival for the test data, 
# and then use it
# as the predictor in a Cox model . We see that the 1st supervised PC is
# highly significant; the last 2 are not
fit.cts <- superpc.predict(train.obj,train.pc,test.pc,threshold=thresh,# What?
                           n.components=3,prediction.type="continuous") 
                          #might want to change threshold and n.components?

temp <- superpc.fit.to.outcome(train.obj, test.pc, fit.cts$v.pred)
#sink("clustCoef.tex")
#  xtable(temp$coef)
#sink()
# Likelihood ratio test=3.76  on 3 df, p=0.288  n= 35, number of events= 20


# sometimes a discrete (categorical) predictor is attractive.
# Here we form two groups by cutting the predictor at its median
# and then plot Kaplan-Meier curves for the two groups
#fit.groups<- superpc.predict(train.obj,train.pc,test.pc,threshold=1,# What?
#                             n.components=1,prediction.type="discrete") 
#                             #might want to change threshold and n.components?
#superpc.fit.to.outcome(train.obj, test.pc, fit.groups$v.pred)


library(rms)
# Doesn't work...
genes <- fit.cts$v.pred[1:nrow(test.exprss),]
colnames(genes) <- c("PC1","PC2","PC3")
surv.pca <- survfit(Surv(test.pc$y, test.pc$censoring.status)~genes[,1]+genes[,2]+genes[,3])
#plot(surv.pca,col=2:3,xlab="time",ylab="Probability of Survival")

# Finally, we look for a predictor of survival a small number of
# markers (rather than all p markers). We do this by computing an importance
# score (see JASA paper, section 4) for each marker
# (if data standardized, equal to its correlation with the first supervised PC predictor).
# Features with the highest correlation contribute
# most to the prediction of the outcome.
# We soft threshold the importance scores, and use the shrunken
# scores as marker weights to form a reduced predictor. Cross-validation
# gives us an estimate of the best amount to shrink and an idea of
# how well the shrunken predictor works.
fit.red<- superpc.predict.red(train.obj, train.pc, test.pc, threshold=thresh)
fit.redcv<- superpc.predict.red.cv(fit.red, cv.obj, train.pc, threshold=thresh)
#pdf("latex/raw/superpcRed.pdf")
#  superpc.plotred.lrtest(fit.redcv)
#dev.off()
#par(mfrow=c(1,1))

# Finally we list the significant markers, in order of decreasing importance score
# raw-score is marginal Cox score
good <- superpc.listfeatures(test.pc, train.obj, fit.red)
good <- good[which(good[,3]!="NA"),]
feat.num <- as.numeric(substr(good[,3],8,nchar(good[,3])))
pca.fit.beta <- rep(0,ncol(train.exprss))

gene <- as.matrix(train.exprss[,feat.num])
pca.res <- coxph(Surv(train.set$time,train.set$cens)~gene)
pca.fit.beta[feat.num] <- pca.res$coef

#pdf("latex/raw/devResPCA.pdf")
#  plot(residuals(pca.res,type="deviance"))
#dev.off()
#
#sink("latex/raw/pcaCoef.tex")
#  xtable(summary(pca.res)$coef)
#  #Likelihood ratio test=59.8  on 18 df, p=2.21e-06  n= 130
#sink()

pca.risk.score.test <- exp(as.matrix(test.exprss)%*%pca.fit.beta)

risk.group <- ifelse(pca.risk.score.test < median(pca.risk.score.test),"Low","High")
#pdf("latex/raw/pcaKM.pdf")
#  survplot(km <- survfit(Surv(test.set$time,test.set$cens)~risk))
#  title(main="PCA KM: High-Risk vs. Low-Risk")
#dev.off()
#lo <- which(risk=="Low")
#hi <- which(risk=="High")
#km.hi <- survfit(Surv(test.set$time[hi],test.set$cens[hi])~risk[hi])
#km.lo <- survfit(Surv(test.set$time[lo],test.set$cens[lo])~risk[lo]) 
#se.pcnt(.5,km.hi)
#se.pcnt(.5,km.lo)
#km

pca.AUC.res <- NA
library(survivalROC)

time.seq <- seq(0,140,by=.5)
for(i in 1:length(time.seq)) {
  pca.AUC.res[i] <- survivalROC(Stime=test.set$time,status=test.set$cens,
                                marker=pca.risk.score.test,
                                predict.time=time.seq[i],method="KM")$AUC
}

#write.table(pca.AUC.res,"aucData/pca.txt",quote=F,col.names=F,row.names=F)


pdf("latex/raw/pcaAUC.pdf")
plot(time.seq,pca.AUC.res,type="l",ylim=c(0,1),main="Time-Dependent ROC",
     ylab="AUC",xlab="Time",col="orange",lwd=3)
dev.off()

