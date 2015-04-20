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





