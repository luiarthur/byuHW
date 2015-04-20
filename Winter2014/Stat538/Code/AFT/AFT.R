
library(survival)

#?survreg
#?survreg.distributions

weib.surv <- function(x,gamma,lambda) {
# lambda = scale = exp(-INTERCEPT)
# gamma = shape = 1/SCALE
  exp(-(x/lambda)^gamma)
}

loglog.surv <- function(x,beta,alpha) {
  1-(((x/alpha)^(-beta)) + 1)^(-1)
}

lognorm.surv <- function(x,mu,sigma) {
  1-pnorm((log(x)-mu)/sigma)
}

# Survival in a randomised trial comparing two treatments for ovarian cancer
# futime:	 survival or censoring time
# fustat:	 censoring status
# age:	         in years
# resid.ds:	 residual disease present (1=no,2=yes)
# rx:	         treatment group
# ecog.ps:	 ECOG performance status (1 is better, see reference)
data(ovarian)

ov.weib <- survreg(Surv(futime,fustat)~1,data=ovarian,dist="weibull")
# lambda = scale = exp(INTERCEPT)
# gamma = shape = 1/SCALE
ov.weib.lambda <- exp(ov.weib$coef[[1]])
ov.weib.gamma <- 1/ov.weib$scale

ov.loglog <- survreg(Surv(futime,fustat)~1,data=ovarian,dist="loglogistic")
# alpha = scale = exp(INTERCEPT)
# beta = shape = 1/SCALE
ov.loglog.alpha <- exp(ov.loglog$coef[[1]])
ov.loglog.beta <- 1/ov.loglog$scale

ov.lognorm <- survreg(Surv(futime,fustat)~1,data=ovarian,dist="lognormal")
ov.lognorm.mu <- ov.lognorm$coef[[1]]
ov.lognorm.sigma <- ov.lognorm$scale


##### nonparametric estimate of survival curve for entire group, compare w/ parametric curve
ov.KM <- survfit(Surv(futime,fustat)~1,data=ovarian,type="kaplan-meier")
summary(ov.KM)

ymax <- 1
ymin <- 0

ov.KM.time <- c(0,ov.KM$time)
ov.KM.surv <- c(1,ov.KM$surv)

# weibull
plot(ov.KM.time,ov.KM.surv,ylim=c(ymin,ymax),type="s",lty=1,
xlab="Time",
ylab="Survival",
main="Ovarian cancer survival",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)
curve(weib.surv(x,ov.weib.gamma,ov.weib.lambda),from=0,to=1200,add=T)

# log-logistic
plot(ov.KM.time,ov.KM.surv,ylim=c(ymin,ymax),type="s",lty=1,
xlab="Time",
ylab="Survival",
main="Ovarian cancer survival",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)
curve(loglog.surv(x,ov.loglog.beta,ov.loglog.alpha),from=0,to=1200,add=T)

# log-normal
plot(ov.KM.time,ov.KM.surv,ylim=c(ymin,ymax),type="s",lty=1,
xlab="Time",
ylab="Survival",
main="Ovarian cancer survival",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)
curve(lognorm.surv(x,ov.lognorm.mu,ov.lognorm.sigma),from=0,to=1200,add=T)


##### fit using covariates
#ov.KM <- survfit(Surv(futime,fustat)~age,data=ovarian,type="kaplan-meier")
ov.loglog <- survreg(Surv(futime,fustat)~age+as.factor(resid.ds)+as.factor(rx)+as.factor(ecog.ps),data=ovarian,dist="loglogistic")
summary(ov.loglog)

### interpretation:
# b_age = -0.0690
# each one-unit increase in age corresponds to a -0.0690 decrease in log survival time
# exp(b_age) = 0.9333
# S_i(t) = S_0(t/exp(b_age))
# - Every additional year in age corresponds to a decrease in survival time by a factor of 0.933
# - Time to death in ovarian cancer patients is accelerated by a factor of 0.933 for every additional year in age
# b_rx = 0.5957
# Survival times for those in treatment group 2 have a predicted 0.5957 increase in log survival time
# exp(b_rx) = 1.814
# S_tx2(t) = S_tx1(t/exp(b_age))
# Survival time is increased by a factor of 1.814 for those assigned to treatment group 2

ov.loglog <- survreg(Surv(futime,fustat)~age,data=ovarian,dist="loglogistic",x=T)
summary(ov.loglog)
alpha <- exp(ov.loglog$coef[[1]])
beta <- 1/ov.loglog$scale



##### plotting

# baseline (reference group)
loglog.surv <- function(x,beta,alpha) {
  1-(((x/alpha)^(-beta)) + 1)^(-1)
}

# incorporate covariates
loglog.surv <- function(t,beta,alpha,xdata,coef) {
  # loglog.surv(exp(coef%*%t(xdata)),beta,alpha)
  1-(((t/(exp(coef%*%t(xdata))*alpha))^(-beta)) + 1)^(-1)
}

xdata <- ov.loglog$x[,-1]
coefmodel <- ov.loglog$coef[-1]

quantile(ovarian$age)

# log-logistic
plot(ov.KM.time,ov.KM.surv,ylim=c(ymin,ymax),type="s",lty=1,
xlab="Time",
ylab="Survival",
main="Ovarian cancer survival",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)
curve(loglog.surv(x,beta,alpha,62,coefmodel),from=0,to=1200,add=T)
curve(loglog.surv(x,beta,alpha,50,coefmodel),from=0,to=1200,add=T)


##### plot separate age groups (divide by 75th percentile)

hist(ovarian$age)
# dichotomize age, rerun model (so will get only two separate age curves)
age.low <- ifelse(ovarian$age < quantile(ovarian$age)[4],1,0)
ov.loglog <- survreg(Surv(futime,fustat)~age.low,data=ovarian,dist="loglogistic")
summary(ov.loglog)

# get parameters corresponding to original formulation of model
alpha <- exp(ov.loglog$coef[[1]])
beta <- 1/ov.loglog$scale
# coefficient(s) corresponding to covariates (exclude intercept, error param)
coefmodel <- ov.loglog$coef[-1]

ymax <- 1
ymin <- 0

ov.agelow.KM <- survfit(Surv(futime[age.low==1],fustat[age.low==1])~1,data=ovarian)
ov.agehigh.KM <- survfit(Surv(futime[age.low==0],fustat[age.low==0])~1,data=ovarian)
summary(ov.agelow.KM)
summary(ov.agehigh.KM)

ov.agelow.KM.time <- c(0,ov.agelow.KM$time)
ov.agelow.KM.surv <- c(1,ov.agelow.KM$surv)
ov.agehigh.KM.time <- c(0,ov.agehigh.KM$time)
ov.agehigh.KM.surv <- c(1,ov.agehigh.KM$surv)

#postscript(file="tumor_KM_compare.ps",height=10,width=10)
plot(ov.agelow.KM.time,ov.agelow.KM.surv,ylim=c(ymin,ymax),type="s",lty=1,
xlab="Time",
ylab="Survival",
main="Ovarian cancer survival",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)
lines(ov.agehigh.KM.time,ov.agehigh.KM.surv,lty=1,type="s",col=2)
curve(loglog.surv(x,beta,alpha,1,coefmodel),from=0,to=1200,add=T)
curve(loglog.surv(x,beta,alpha,0,coefmodel),from=0,to=1200,add=T,col=2)



#####################################################
### MODEL FIT

ov.loglog <- survreg(Surv(futime,fustat)~age,data=ovarian,dist="loglogistic")

#### assessing AFT assumption: percentile-percentile plots
# percentile-percentile plot for age
# (need to do separately for each covariate)

g1.KM <- ov.agelow.KM 
g2.KM <- ov.agehigh.KM

g1.perc <- g2.perc <- NA
p <- seq(.1,1,by=.1)
for(i in 1:length(p)) {
  g2.perc[i] <- min(g2.KM$time[g2.KM$surv <= (1-p[i])])
  g1.perc[i] <- min(g1.KM$time[g1.KM$surv <= (1-p[i])])
}

index <- min(c(sum(g2.perc < Inf), sum(g1.perc < Inf)))
plot(g2.perc[1:index],g1.perc[1:index],type="l")



#### assessing parametric assumptions
# probability plots
summary(ov.KM)


# weibull
plot(log(ov.KM$time), log(log(1/(ov.KM$surv))),pch=19)

# lognormal
plot(log(ov.KM$time), qnorm(1-ov.KM$surv),pch=19)

# loglogistic
plot(log(ov.KM$time), log((1/(ov.KM$surv))-1),pch=19)




#### assessing overall fit: deviance residuals
ov.loglog.dev <- residuals(ov.loglog,type="deviance")
plot(ov.loglog.dev,pch=19)


