
# Survival in patients with lung cancer at Mayo Clinic. Performance scores rate
# how well the patient can perform usual daily activities.

#       inst:       Institution code
#       time:       Survival time in days
#       status:     censoring status 1=censored, 2=dead
#       age:        Age in years
#       sex:        Male=1 Female=2
#       ph.ecog:    ECOG performance score (0=good 5=dead)
#       ph.karno:   Karnofsky performance score (bad=0-good=100) rated by physician
#       pat.karno:  Karnofsky performance score  rated by patient
#       meal.cal:   Calories consumed at meals
#       wt.loss:    Weight loss in last six months

library(survival)
lungData <- lung
lungData$status <- ifelse(lungData$status==2,1,0)

set.seed(3)
trainIndex <- sample(c(1:228),114)
lungDataTrain <- lungData[trainIndex,]
lungDataTest <- lungData[-trainIndex,]
head(lungDataTrain)



# Model 1
coxph(Surv(time,status) ~ as.factor(sex) + ph.ecog, data=lungData,x=T)
lungCoxMod1 <- coxph(Surv(time,status) ~ as.factor(sex) + ph.ecog ,data=lungDataTrain,x=T)
summary(lungCoxMod1)
lungCoxCoefMod1 <- c(0,as.numeric(lungCoxMod1$coef),0,0,0,0)


# Model 2
coxph(Surv(time,status) ~ age, data=lungData,x=T)
lungCoxMod2 <- coxph(Surv(time,status) ~ age,data=lungDataTrain,x=T)
summary(lungCoxMod2) 
lungCoxCoefMod2 <- c(as.numeric(lungCoxMod2$coef),rep(0,6))



###### ASSESS PERFORMANCE OF TWO APPROACHES: RISK GROUPS ##########################

## risk-groups (check ability of two approaches to more accurately identify
##   high- and low-risk groups
##   make two-group comparisions, check via plots and log-rank test)


# risk score for those in training data set
exp(lungCoxCoefMod1%*%t(as.matrix(lungDataTrain[,c(4:10)])))
# risk score for those in test data set
exp(lungCoxCoefMod1%*%t(as.matrix(lungDataTest[,c(4:10)])))
riskScoreTestMod1 <- exp(lungCoxCoefMod1%*%t(as.matrix(lungDataTest[,c(4:10)])))
riskScoreTestMod2 <- exp(lungCoxCoefMod2%*%t(as.matrix(lungDataTest[,c(4:10)])))

# create test set risk groups based on risk score (predicted outcome)
# 0: low-risk group
# 1: high-risk group
riskGroupMod1 <- ifelse(riskScoreTestMod1 < median(riskScoreTestMod1,na.rm=T),0,1)
riskGroupMod2 <- ifelse(riskScoreTestMod2 < median(riskScoreTestMod2,na.rm=T),0,1)

# check whether penalized approach identified markers that
# successfully distinguish between high-risk / low-risk groups
survdiff(Surv(time,status)~c(riskGroupMod1),data=lungDataTrain)
survdiff(Surv(time,status)~c(riskGroupMod2),data=lungDataTrain)

survdiff(Surv(time,status)~c(riskGroupMod1),data=lungDataTest)
survdiff(Surv(time,status)~c(riskGroupMod2),data=lungDataTest)


# plot K-M survival functions for high/low risk (full group)

par(mfrow=c(1,2))
g1.time <- lungDataTest$time[riskGroupMod1==0] # low-risk group
g1.cens <- lungDataTest$status[riskGroupMod1==0]

g2.time <- lungDataTest$time[riskGroupMod1==1] # high-risk group
g2.cens <- lungDataTest$status[riskGroupMod1==1]

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

plot(g1.KM.time,g1.KM.surv,ylim=c(ymin,ymax),type="s",lty=1,
     xlab="Time",
     ylab="Survival",
     main="Test: Model 1",
     cex.main=2,
     cex.axis=1.5,
     cex.lab=1.5,
     xlim=c(0,max(lungDataTest$time)))
lines(g2.KM.time,g2.KM.surv,lty=1,type="s",col=2)




g1.time <- lungDataTest$time[riskGroupMod2==0] # low-risk group
g1.cens <- lungDataTest$status[riskGroupMod2==0]

g2.time <- lungDataTest$time[riskGroupMod2==1] # high-risk group
g2.cens <- lungDataTest$status[riskGroupMod2==1]

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

plot(g1.KM.time,g1.KM.surv,ylim=c(ymin,ymax),type="s",lty=1,
     xlab="Time",
     ylab="Survival",
     main="Test: Model 2",
     cex.main=2,
     cex.axis=1.5,
     cex.lab=1.5,
     xlim=c(0,max(lungDataTest$time)))
lines(g2.KM.time,g2.KM.surv,lty=1,type="s",col=2)




###### ASSESS PERFORMANCE OF TWO APPROACHES: TIME-DEPENDENT ROC ##########################

library(survivalROC)

lungDataTest.time <- lungDataTest$time
lungDataTest.status <- lungDataTest$status

time.seq <- seq(0,1010,by=1)


marginal.AUC.res <- NA
for(i in 1:length(time.seq)) {
marginal.AUC.res[i] <- survivalROC(Stime=lungDataTest.time,status=lungDataTest.status,marker=riskScoreTestMod1,predict.time=time.seq[i],method="KM")$AUC
}
marginal.AUC.res

pen.AUC.res <- NA
for(i in 1:length(time.seq)) {
pen.AUC.res[i] <- survivalROC(Stime=lungDataTest.time,status=lungDataTest.status,marker=riskScoreTestMod2,predict.time=time.seq[i],method="KM")$AUC
}
pen.AUC.res

par(mfrow=c(1,1))
plot(marginal.AUC.res,type="l",
     ylim=c(0,1),
     main="Time-Dependent ROC",
     ylab="AUC",
     xlab="Time")
lines(pen.AUC.res,type="l",col="red")




