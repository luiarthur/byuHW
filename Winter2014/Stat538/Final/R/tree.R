library(survival)
library(randomForestSRC) # updated version of randomSurivalForest
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


  library(foreach)
  library(doMC)
  registerDoMC(16)

  one.summary <- function(i) {
    train.summary <- summary(coxph(Surv(train.set$time,train.set$cens) ~ 
                             train.exprss[,i]))
    p <- train.summary$coef[5]
    z <- train.summary$coef[4]
    matrix(c("pval"=p,"zscore"=z),1,2)
  }

  system.time(
    pz <- foreach(i=1:(ncol(train.exprss)),.combine=rbind) %dopar% 
                  one.summary(i) # Takes 35s
  )  

  colnames(pz) <- c("pval","zscore")
  pz <- data.frame(pz)

  pvalue <- pz[,1]
  zscore <- pz[,2]
  

sum(pvalue < 1e-5)
(marker.index <- which(pvalue < 1e-5))
subset <- train.exprss[,marker.index]


# ntree: number of bootstrap samples to be drawn from original data set
# mtry: number of variables randomly sampled at each split - default is sqrt(p)
# nodesize: minimum number of deaths with unique survival times required for a terminal node. Default is
#      approximately 3 for right-censored data
# splitrule: Four primary splitting rules are available for growing a survival forest for right-censored data:
#      "logrank", "conserve", "logrankscore", and "random".
#      The default rule, "logrank", splits tree nodes by maximization of the log-rank test statistic.

df <- as.data.frame(subset)
df$time <- train.set$time
df$time <- ifelse(df$time==0,.01,df$time)
df$cens <- train.set$cens
v.out <- rfsrc(Surv(time,cens) ~ .,ntree=500,forest=T,data=df,mtry=10)
#print(v.out)



#pdf("latex/raw/varImp.pdf",height=11)
#  plot(v.out)
#dev.off()
#system("firefox latex/raw/tree.varImp.pdf")

# The error rate is smaller than 0.5, the
# benchmark value associated with a procedure no better
# than flipping a coin. This is evidence
# that Karnofsky score is predictive.
# We can investigate the effect of Karnofsky score
# more closely by considering how the ensemble estimated
# mortality varies as a function of the predictor:

# takes awhile
#postscript("rfs.ps",height=10,width=10)
#pdf("latex/raw/tree.pdf")
#  plot.variable(v.out, partial=T) # Takes 10 minutes
#dev.off()
#system("ps2pdf rfs.ps rfs.pdf")
#system("rm rfs.ps")

# Vertical axis is mortality for a given Karnofsky value x and
# represents expected number of deaths.


##########################################################################
##### construct Cox model using identified variables with high importance

v.out$importance
quantile(v.out$importance)
sum(v.out$importance > .002)
marker.index[v.out$importance > .002]
final.marker.index <- marker.index[v.out$importance > .002]

genes <- as.matrix(train.exprss[,final.marker.index])
rfs.cox <- coxph(Surv(train.set$time,train.set$cens)~genes)
library(xtable)
#sink("latex/raw/treeCox.tex")
#  xtable(rfs.cox)
#sink()  
summary(rfs.cox)

# need to create a coefficient vector (length=number of markers,
#    all zeros expect for selected markers
rfs.fit.beta <- rep(0,ncol(train.exprss))
for(i in 1:length(final.marker.index)) {
  rfs.fit.beta[final.marker.index[i]] <- rfs.cox$coef[[i]]
}

# risk score for those in test data set
rfs.risk.score.test <- exp(as.matrix(test.exprss)%*%rfs.fit.beta)

# create test set risk groups based on risk score (predicted outcome)
# 0: low-risk group
# 1: high-risk group
risk.group <- ifelse(rfs.risk.score.test < median(rfs.risk.score.test),0,1)

# check whether penalized approach identified markers that
# successfully distinguish between high-risk / low-risk groups
# plot K-M survival functions for high/low risk

library(rms)
risk <- ifelse(risk.group==1,"High","Low")
#pdf("latex/raw/treeKM.pdf")
#  survplot(km <- survfit(Surv(test.set$time,test.set$cens)~risk))
#  title(main="Random Forest: High-Risk vs. Low-Risk")
#dev.off()

hi <- which(risk=="High")
km.hi <- survfit(Surv(test.set$time[hi],test.set$cens[hi])~risk[hi,])

lo <- which(risk=="Low")
km.lo <- survfit(Surv(test.set$time[lo],test.set$cens[lo])~risk[lo,])

cox.hilo <- coxph(Surv(test.set$time,test.set$cens)~risk)
# logrank test: .02, 1 df => p = 0.8958
#pdf("latex/raw/rfResid.pdf")
#  plot(residuals(rfs.cox,type="deviance"),pch=20,col="green")
#dev.off()

#sink("latex/raw/treeKM.tex")
#  tab <- rbind(se.pcnt(.5,km.lo)$CI,se.pcnt(.5,km.hi)$CI)
#  rownames(tab) <- c("Low-Risk Group","High-Risk Group")
#  xt <- xtable(tab)
#  caption(xt) <- "\\tiny Median"
#  print(xt,caption.placement="top")
#sink()
#km

library(survivalROC)

time.seq <- seq(0,140,by=.5)

rf.AUC.res <- NA
for(i in 1:length(time.seq)) {
  rf.AUC.res[i] <- survivalROC(Stime=test.set$time,status=test.set$cens,
                               marker=rfs.risk.score.test,predict.time=time.seq[i],
                               method="KM")$AUC
}
#write.table(rf.AUC.res,"aucData/rf.txt",quote=F,col.names=F,row.names=F)

#pdf("latex/raw/rfAUC.pdf")
#  plot(time.seq,rf.AUC.res,type="l",ylim=c(0,1),main="Time-Dependent ROC",
#       ylab="AUC",xlab="Time",col="green",lwd=3)
#dev.off()


