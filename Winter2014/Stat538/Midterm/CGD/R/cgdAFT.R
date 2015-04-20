rm(list=ls())
library(survival)

# center
# treat
# sex
# age
# height
# weight
# inherit
# steroids
# propylac
# hos.cat
# time (days)
# status (1 = infection, 0 = censored) Want dead=infection=1 => data is good
cgd <- read.csv("../../Data/cgd.csv")[,-c(1:2)]

#0: survival functions
sweib <- function(x,gamma,lambda) exp(-(x/lambda)^gamma)
sllog <- function(x,beta,alpha) 1-(((x/alpha)^(-beta)) + 1)^(-1)
slnorm <- function(x,mu,sigma) pnorm((log(x)-mu)/sigma,lower=F)


#1: Fit AFT Models: (Weibull, loglogistic, lognormal)
weib <- survreg(Surv(time,status)~1,data=cgd,dist="weibull")
weib.lambda <- exp(weib$coef[[1]]) # lambda = scale = exp(Intercept)
weib.gamma <- 1/weib$scale         # gamma = shape = 1/scale

llog <- survreg(Surv(time,status)~1,data=cgd,dist="loglogistic")
llog.lambda <- exp(llog$coef[[1]]) # lambda = scale = exp(Intercept)
llog.gamma <- 1/llog$scale         # gamma = shape = 1/scale

lnorm <- survreg(Surv(time,status)~1,data=cgd,dist="lognormal")
lnorm.mu <- lnorm$coef[[1]]
lnorm.sigma <- lnorm$scale

#2: Plot Parametric and Non-parametric (KM) Curves:
KM <- survfit(Surv(time,status)~1,data=cgd,type="kaplan-meier")

plot.KM <- function(cex=1,compare=T){
  plot(c(0,KM$time),c(1,KM$surv),ylim=c(0,1),type="s",lty=1,
       xlab="Time to Infection",ylab="Survival",
       main="Survival Curve for Time to Infection",lwd=3)
  if (compare) {
  to <- max(KM$time)
    curve(sweib(x,weib.gamma,weib.lambda),fr=0,to=to,add=T,col="blue",lwd=2)
    curve(sllog(x,llog.gamma,llog.lambda),fr=0,to=to,add=T,col="red",lwd=2)
    curve(slnorm(x,lnorm.mu,lnorm.sigma),fr=0,to=to,add=T,col="green",lwd=2)
    legend("topright",legend=c("K.M. Curve","Weibull","loglogistic","lognormal"),
           col=c("black","blue","red","green"),lwd=3,cex=cex)
  # It appears that the Weibull curve follows the K.M. curve the closest.       
  }
}


#3: Assess Parametric Assumptiosn:
    #Probability plots
    summary(KM)

    # Weibull
    plot.weib.prob <- function()
      plot(log(KM$time), log(log(1/(KM$surv))),pch=19,
           xlab="Log of Time to Infection",
           ylab="log of Hazard of Time to Infection",
           main="Probability Plot for Weibull",col="blue")

    # Lognormal
    plot.lnorm.prob <- function()
      plot(log(KM$time), qnorm(1-KM$surv),pch=19,
           xlab="Log of Time to Infection",
           ylab="Percentiles of Time to Infection",
           main="Probability Plot for Lognormal",col="green")

    # Loglogistic
    plot.llog.prob <- function()
      plot(log(KM$time), log((1/(KM$surv))-1),pch=19,
           xlab="Log of Time to Infection",
           ylab="Negative Log Odds of Time to Infection",
           main="Probability Plot for Loglogistic",col="red")
    
    # They all look bad. But go with Weibull because plot.KM appears to be closest.


#4: Fit using Covariates:

  # Parametric Model:

  age.cut <- quantile(cgd$age,.42)
  cgd.new <- cgd
  cgd.new$age <- ifelse(cgd$age>=age.cut,1,0)
  cgd.new$treat <- ifelse(cgd$treat=="placebo",0,1)
  weib.full <- survreg(Surv(time,status)~.,dat=cgd.new,dist="weibull")
  weib.1    <- survreg(Surv(time,status)~1,dat=cgd.new,dist="weibull")
  weib.step <- step(weib.1,scope=list(lower=weib.1,upper=weib.full),
                    data=cgd.new,direction="both")

  # Parameter Estimates:
  l <- exp(weib.step$coef[[1]]) # lambda = scale = exp(Intercept)
  g <- 1/weib.step$scale         # gamma = shape = 1/scale
  coef <- weib.step$coef[-1]

  # Summary Table:
  summary(weib.step)

  KM.treat <- survfit(Surv(time,status)~treat,data=cgd,type="kaplan-meier")

  plot.tmt.s <- function(cex=1,age.cutoff=age.cut,lwd=1) {
    plot(KM.treat,xlab="Time",ylab="Survival",
         main="Survival Curve for Time to Infection",
         lty=1,lwd=lwd,col=c("purple","red")) #The first one is the 0

    to <- max(KM.treat$time)
    eeta.0.0 <- exp(c(0,0) %*% coef) # treatment (placebo), age (<10)
    eeta.1.0 <- exp(c(1,0) %*% coef)
    eeta.0.1 <- exp(c(0,1) %*% coef)
    eeta.1.1 <- exp(c(1,1) %*% coef) 

    curve(sweib(x/eeta.0.0,g,l),fr=0,to=to,add=T,col='purple',lwd=lwd,lty=3)
    curve(sweib(x/eeta.1.0,g,l),fr=0,to=to,add=T,col='red',   lwd=lwd,lty=3)
    curve(sweib(x/eeta.0.1,g,l),fr=0,to=to,add=T,col='purple',lwd=lwd,lty=2)
    curve(sweib(x/eeta.1.1,g,l),fr=0,to=to,add=T,col='red',   lwd=lwd,lty=2)

    legend("bottomleft",col=c("Purple","Red"),lty=c(1,1,2,2,3,3),
           leg=c(paste(c("Placebo","rIFN-g"),"KM Survival Surve"),
                 paste(c("Placebo","rIFN-g"),"Survival Curve at age >=",age.cutoff),
                 paste(c("Placebo","rIFN-g"),"Survival Curve at age <",age.cutoff)),
                 cex=cex,lwd=lwd)
    # People age 10 or older with treatment survive longer
  }
  
#5: Model Fit:

  # Assess AFT Assumption: Q-Q plot
  # One for each covariate

  cgd.ty <- cgd.new[which(cgd.new$treat==1 & cgd.new$age==0),]
  cgd.to <- cgd.new[which(cgd.new$treat==1 & cgd.new$age==1),]
  cgd.py <- cgd.new[which(cgd.new$treat==0 & cgd.new$age==0),]
  cgd.po <- cgd.new[which(cgd.new$treat==0 & cgd.new$age==1),]

  KM.ty <- survfit(Surv(time,status)~1,data=cgd.ty)
  KM.to <- survfit(Surv(time,status)~1,data=cgd.to)
  KM.py <- survfit(Surv(time,status)~1,data=cgd.py)
  KM.po <- survfit(Surv(time,status)~1,data=cgd.po)

  ty.p <- to.p <- py.p <- po.p <- NA
  #p <- seq(.1,1,by=.1)
  p <- 1:100/100
  for(i in 1:length(p)) {
    ty.p[i] <- min(KM.ty$time[KM.ty$surv <= (1-p[i])])
    to.p[i] <- min(KM.to$time[KM.to$surv <= (1-p[i])])
    py.p[i] <- min(KM.py$time[KM.py$surv <= (1-p[i])])
    po.p[i] <- min(KM.po$time[KM.po$surv <= (1-p[i])])
  }

  t.index <- min(c(sum(ty.p < Inf), sum(ty.p < Inf)))
  plot.t.qq <- function(){
    plot(ty.p[1:t.index],ty.p[1:t.index],type="p",
         main="Q-Q Plot for rIFN-g",col="red",lwd=3,
         xlab="Percentile for rIFN-g age < 10 Group",
         ylab="Percentile for rIFN-g age >= 10 Group")
  }

  p.index <- min(c(sum(py.p < Inf), sum(po.p < Inf)))
  plot.p.qq <- function(){
    plot(py.p[1:p.index],po.p[1:p.index],type="p",
         main="Q-Q Plot for Females",col="red",lwd=3,
         xlab="Percentile for Placebo age < 10 Group",
         ylab="Percentile for Placebo age >= 10 Group")
  }

#6: Assessing Overall Fit: Deviance Residuals
    weib.dev <- residuals(weib.step,type="deviance")
    plot.resid <- function()
      plot(weib.dev,pch=19,col='purple',main="Residuals Plot")

# plots: ###########################################################
  plot.KM()
  plot.tmt.s(cex=.9,lwd=1.3)
  plot.weib.prob()
  plot.lnorm.prob()
  plot.llog.prob()
  plot.t.qq()
  plot.p.qq()
  plot.resid()

  table(cgd$treat,cgd$status) # placebo has no infect on infection, tmt does.
  table(cgd.new$age,cgd$status) # high age reduces infection. 

