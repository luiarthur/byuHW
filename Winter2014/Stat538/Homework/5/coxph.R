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

data.lung$status <- ifelse(data.lung$status==2,1,0)

data.lung <- lung
head(data.lung)
str(data.lung)
hist(data.lung$age)
hist(data.lung$sex)
hist(data.lung$ph.ecog)
hist(data.lung$ph.karno)
hist(data.lung$meal.cal)
hist(data.lung$wt.loss)
  
data.lung$sex <- data.lung$sex-1

#lung1 <- lung[1:114,]
#lung2 <- lung[155:228,]

# full model
lung.res <- coxph(Surv(time,status)~age + as.factor(sex) + ph.ecog + ph.karno + pat.karno + meal.cal + wt.loss,data=data.lung,x=T)
summary(lung.res)

# reduced model
lung.res.reduced <- coxph(Surv(time,status)~ sex + ph.ecog + ph.karno + pat.karno + wt.loss,data=data.lung,x=T)
summary(lung.res.reduced)

# gender model
lung.res.sex <- coxph(Surv(time,status)~ sex,data=data.lung,x=T)
summary(lung.res.sex)

# null model
lung.res.null <- coxph(Surv(time,status)~1,data=data.lung,x=T)
summary(lung.res.null)

# bad model (for residual analysis)
lung.res.bad <- coxph(Surv(time,status)~meal.cal + wt.loss,data=data.lung,x=T)
summary(lung.res.bad)



# get loglik for cox ph model
# likelihood under alternative (reduced) / likelihood under null (full)
#      reject alternative for small values of this ratio
#      (for large negative values of log ratio) 
#      maximum is optimal
#      first value is that obtained using initial values (ignore)
#      second value is that obtained using final coefficient values (use)
lung.res$loglik
lung.res.reduced$loglik
lung.res.null$loglik

#### Likelihood Ratio Test
# log likelihood for reduced model
test.stat <- -2*(lung.res.reduced$loglik[2] - lung.res$loglik[2])
test.stat
1-pchisq(test.stat,df=7-5)
# reduced model is rejected

test.stat <- -2*(lung.res.null$loglik[1] - lung.res.reduced$loglik[2])
test.stat
1-pchisq(test.stat,df=5)
# null model is rejected (in comparison to reduced)

test.stat <- -2*(lung.res.null$loglik[1] - lung.res$loglik[2])
test.stat
1-pchisq(test.stat,df=7)
# null model is rejected (in comparison to full)



#### get AIC for cox ph model
# AIC = -2loglike + 2p, where p is number of parameters in model
# minimum is optimal
extractAIC(lung.res)
extractAIC(lung.res.reduced)



# ASSUMPTION OF LIKELIHOOD-BASED METHODS: data sets are the same
# this is a problem






########## Plotting Results

summary(lung.res.reduced)
# want to look at gender, ph.ecog groups
hist(lung$ph.ecog)
ph.ecog.dichot <- ifelse(lung$ph.ecog <= 1,0,1)

g1.KM <- survfit(Surv(time[ph.ecog.dichot==0 & sex==1],status[ph.ecog.dichot==0 & sex==1])~1,type="kaplan-meier",data=lung)
g2.KM <- survfit(Surv(time[ph.ecog.dichot==0 & sex==2],status[ph.ecog.dichot==0 & sex==2])~1,type="kaplan-meier",data=lung)
g3.KM <- survfit(Surv(time[ph.ecog.dichot==1 & sex==1],status[ph.ecog.dichot==1 & sex==1])~1,type="kaplan-meier",data=lung)
g4.KM <- survfit(Surv(time[ph.ecog.dichot==1 & sex==2],status[ph.ecog.dichot==1 & sex==2])~1,type="kaplan-meier",data=lung)

summary(g1.KM)
summary(g2.KM)
summary(g3.KM)
summary(g4.KM)

ymax <- 1
ymin <- 0

g1.KM.time <- c(0,g1.KM$time)
g2.KM.time <- c(0,g2.KM$time)
g3.KM.time <- c(0,g3.KM$time)
g4.KM.time <- c(0,g4.KM$time)

g1.KM.surv <- c(1,g1.KM$surv)
g2.KM.surv <- c(1,g2.KM$surv)
g3.KM.surv <- c(1,g3.KM$surv)
g4.KM.surv <- c(1,g4.KM$surv)

max(g1.KM.time)
max(g2.KM.time)
max(g3.KM.time)
max(g4.KM.time)


#postscript(file="lung.sex.phecog.ps",height=10,width=10)
plot(g3.KM.time,g3.KM.surv,ylim=c(ymin,ymax),xlim=c(0,1000),type="s",lty=1,
xlab="Time (days)",
ylab="Survival",
main="Lung Cancer Survival",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5,
     )
lines(g1.KM.time,g1.KM.surv,lty=1,type="s",col="red")
lines(g2.KM.time,g2.KM.surv,lty=1,type="s",col="blue")
lines(g4.KM.time,g4.KM.surv,lty=1,type="s",col="green")
legend(580,1.02,"low ph.ecog, male",lty=1,col="red",bty="n")
legend(580,.98,"high ph.ecog, male",lty=1,col="black",bty="n")
legend(580,.94,"low ph.ecog, female",lty=1,col="blue",bty="n")
legend(580,.90,"high ph.ecog, female",lty=1,col="green",bty="n")
#dev.off()



########## Estimation of baseline hazard and baseline survival

survfit(lung.res.reduced)
summary(survfit(lung.res.reduced)) # baseline hazard
coxph.detail(lung.res.reduced)$hazard # hazard increment

# h_i(t) = h_0(t)exp(x'B)
# H_i(t) = int_0^t h_i(t)
#        = int_0^t h_0(t)exp(x'B)
#        = exp(x'B) int_0^t h_0(t)
#        = exp(x'B) H_o(t)
# H_0(t) = -log(S_o(t))
# H_i(t) = -log(S_i(t))
# S_o(t) = exp(-H_o(t))
# S_i(t) = exp(-H_i(t))
#        = exp(-exp(x'B)H_o(t))
#        = exp(-H_o(t))^exp(x'B)
#        = S_o(t)^exp(x'B)
# log(S_i(t)) = exp(x'B)log(S_o(t))
# log(S_i(t))/log(S_o(t)) = exp(x'B)






############################################



########## Martingale residuals
lung.res$residuals
residuals(lung.res,type="martingale") # default

plot(c(1:length(lung.res$residuals)),lung.res$residuals,xlab="Index",ylab="Martingale Residuals")

########## Deviance residuals
residuals(lung.res,type="deviance")

plot(c(1:length(residuals(lung.res,type="deviance"))),residuals(lung.res,type="deviance"),xlab="Index",ylab="Deviance Residuals")




########## Deviance residuals
residuals(lung.res,type="deviance")

par(mfrow=c(1,2))

# full model
plot(c(1:length(residuals(lung.res,type="deviance"))),residuals(lung.res,type="deviance"),xlab="Index",ylab="Deviance Residuals")

# reduced model
plot(c(1:length(residuals(lung.res.reduced,type="deviance"))),residuals(lung.res.reduced,type="deviance"),xlab="Index",ylab="Deviance Residuals")

##### to plot against various indices (e.g., time, covariate values),
#####       need to identify those individuals excluded from the model
#####       and create vector consisting of variables of interest (minus those
#####       excluded from the model)

### plot against time

par(mfrow=c(1,1))

# reduced data consist of a limited number of covariates 
data.lung.reduced <- data.lung[,c(2,3,5:8,10)]

# eliminate individuals with missing covariate data
# (R already did this in building model, don't have residuals for these individuals)
data.lung.nona.time <- data.lung.reduced$time[apply(apply(data.lung.reduced,2,is.na),1,sum)==0]

plot(data.lung.nona.time,residuals(lung.res.reduced,type="deviance"),xlab="Time",ylab="Deviance Residuals")

### plot against specific covariate values

data.lung.nona.sex <- data.lung.reduced$sex[apply(apply(data.lung.reduced,2,is.na),1,sum)==0]

plot(data.lung.nona.sex,residuals(lung.res.reduced,type="deviance"),xlab="Sex",ylab="Deviance Residuals")

### plot against specific covariate values

data.lung.nona.phkarno <- data.lung.reduced$ph.karno[apply(apply(data.lung.reduced,2,is.na),1,sum)==0]

plot(data.lung.nona.phkarno,residuals(lung.res.reduced,type="deviance"),xlab="Age",ylab="Deviance Residuals")

### plot against risk scores

score <- lung.res.reduced$coef%*%t(lung.res.reduced$x)
plot(score,residuals(lung.res.reduced,type="deviance"),xlab="Score",ylab="Deviance Residuals")



##### Evidence of poor fit: poor covariate selection (lung cancer data)
########## Martingale residuals
plot(c(1:length(lung.res.bad$residuals)),lung.res.bad$residuals,xlab="Index",ylab="Martingale Residuals",pch=19,cex=.5)
### plot against time
data.lung.reduced <- data.lung[,c(2,3,9,10)]
data.lung.nona.time <- data.lung.reduced$time[apply(apply(data.lung.reduced,2,is.na),1,sum)==0]
plot(data.lung.nona.time,lung.res.bad$residuals,xlab="Time",ylab="Martingale Residuals",pch=19,cex=.5)

########## Deviance residuals
plot(c(1:length(residuals(lung.res.bad,type="deviance"))),residuals(lung.res.bad,type="deviance"),xlab="Index",ylab="Deviance Residuals",pch=19,cex=.5)
### plot against time
plot(data.lung.nona.time,residuals(lung.res.bad,type="deviance"),xlab="Time",ylab="Deviance Residuals",pch=19,cex=.5)

########## test proportional hazards assumption
# significance indicates non-proportionality
cox.zph(lung.res.bad)
# plot scaled Schoenfeld residuals along a smooth curve
postscript("lung.res.bad.ps",height=10,width=10)
plot(cox.zph(lung.res.bad))
dev.off()
system("ps2pdf lung.res.bad.ps lung.res.bad.pdf")
system("rm lung.res.bad.ps")




########## Proportional hazards assumption 
#curve(dweibull(x,5,2),from=0,to=10)
dat1 <- rweibull(200,20,60)
dat2 <- rweibull(200,1,60)
dat.time <- c(dat1,dat2)
dat.cens <- rbinom(length(dat.time),1,.8)
dat.tx <- c(rep(0,length(dat1)),rep(1,length(dat2)))

sim.coxph.plot <- coxph(Surv(dat.time,dat.cens)~strata(dat.tx))
sim.coxph <- coxph(Surv(dat.time,dat.cens)~dat.tx)
summary(sim.coxph)

### plotting the survival functions #######
plot(survfit(sim.coxph.plot))

########## test proportional hazards assumption for simulated data
residuals(lung.res,type="schoenfeld")

# plot scaled Schoenfeld residuals along a smooth curve
plot(cox.zph(sim.coxph))
# significance indicates non-proportionality
cox.zph(sim.coxph)



########## test proportional hazards assumption for lung data (full model)
# plot scaled Schoenfeld residuals along a smooth curve
postscript("lung.res.ps",height=10,width=10)
plot(cox.zph(lung.res))
dev.off()
system("ps2pdf lung.res.ps lung.res.pdf")
system("rm lung.res.ps")

# plot scaled Schoenfeld residuals along a smooth curve
postscript("lung.res.reduced.ps",height=10,width=10)
plot(cox.zph(lung.res.reduced))
dev.off()
system("ps2pdf lung.res.reduced.ps lung.res.reduced.pdf")
system("rm lung.res.reduced.ps")

# significance indicates non-proportionality
cox.zph(lung.res)
cox.zph(lung.res.reduced)








########## Stratified Cox model (use when proportional hazards assumption not met)
# Stratified Cox model provides separate baseline hazard function for each stratum level.
# When the effect of a covariate is found to be time-varying (through use of Schoenfeld residuals),
# can use stratified Cox model

### stratify on meal.cal (had problem in full model)

# NAs present, need to use na.rm
lung$meal.cal <= median(lung$meal.cal)
lung$meal.cal <= median(lung$meal.cal,na.rm=T)
mealcal.dichot <- ifelse(lung$meal.cal <= median(lung$meal.cal,na.rm=T),0,1)

lung.res.st <- coxph(Surv(time,status) ~ strata(mealcal.dichot) + age + as.factor(sex) + ph.ecog + ph.karno + pat.karno + wt.loss,data=data.lung,x=T)
summary(lung.res.st)

cox.zph(lung.res.st)

# plotting the survival functions #######

plot(survfit(lung.res.st))
# doesn't plot different colors, line types for genders

summary(survfit(lung.res.st))$strata
strata.res <- summary(survfit(lung.res.st))$strata
n.mealcal0 <- sum(strata.res=="mealcal.dichot=0")
n.mealcal1 <- sum(strata.res=="mealcal.dichot=1")

mealcal0.time <- c(0,summary(survfit(lung.res.st))$time[1:n.mealcal0])
mealcal1.time <- c(0,summary(survfit(lung.res.st))$time[-c(1:n.mealcal0)])

mealcal0.surv <- c(1,summary(survfit(lung.res.st))$surv[1:n.mealcal0])
mealcal1.surv <- c(1,summary(survfit(lung.res.st))$surv[-c(1:n.mealcal0)])

ymax <- 1
ymin <- 0

#postscript(file="lung.mealcalstrata.ps",height=10,width=10)
plot(mealcal0.time,mealcal0.surv,ylim=c(ymin,ymax),type="s",lty=2,
  xlab="Time", ylab="Survival",
  main="Lung Cancer Survival Curves (Low vs. High Calories per Meal)",
  cex.main=2, cex.axis=1.5, cex.lab=1.5)

lines(mealcal1.time,mealcal1.surv,lty=1,type="s")
legend(580,1.02,"low meal cal",lty=2,bty="n")
legend(580,.98,"high meal cal",lty=1,bty="n")
#dev.off()





########## check for influential observations

residuals(lung.res.reduced,type="score")
dim(residuals(lung.res.reduced,type="score"))
lung.res.reduced$var

getlikdis <- function(x) {
  likdis <- NA
  for (i in 1:nrow(residuals(x,type="score"))) {
    likdis[i] <- t(as.vector(residuals(x,type="score")[i,])) %*% x$var %*%
      as.vector(residuals(x,type="score")[i,])
  }
  return(likdis)
}
as.vector(residuals(lung.res.reduced,type="score")[1,])
t(as.vector(residuals(lung.res.reduced,type="score")[1,])) %*% lung.res.reduced$var %*%
  as.vector(residuals(lung.res.reduced,type="score")[1,])

  
### obtain for candidate model

likdis.lung.res <- getlikdis(lung.res.reduced) #need this


data.lung.reduced <- data.lung[,c(2,3,5:8,10)]
data.lung.nona.time <- data.lung.reduced$time[apply(apply(data.lung.reduced,2,is.na),1,sum)==0]
data.lung.reduced.nona <- data.lung.reduced[apply(apply(data.lung.reduced,2,is.na),1,sum)==0,]


plot(data.lung.nona.time,likdis.lung.res,xlab="Time",ylab="Likelihood Displacement") # need this


### eliminate those with high LD values, check effect on candidate model

data.lung.new <- as.data.frame(data.lung.reduced.nona[likdis.lung.res < .15,])

lung.res.reduced.new <- coxph(Surv(time,status)~ sex + ph.ecog + ph.karno + pat.karno + wt.loss,data=data.lung.new,x=T)
lung.res.reduced.new
lung.res.reduced

likdis.lung.res.new <- getlikdis(lung.res.reduced.new)
data.lung.nona.time.new <- data.lung.new$time[apply(apply(data.lung.new,2,is.na),1,sum)==0]
plot(data.lung.nona.time.new,likdis.lung.res.new,xlab="Time",ylab="Likelihood Displacement")


