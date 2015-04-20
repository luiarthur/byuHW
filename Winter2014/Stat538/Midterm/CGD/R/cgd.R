rm(list=ls())
library(survival)
library(rms) # For: survplot

# My own function :)
se.pcnt <- function(p,km,e=.05){# p is a percentage for a percentile
  # Computes the standard error of a given percentile
  se.S <- km$std.err[which(km$surv<p)[1]]
  perc <-  km$time[which(km$surv<p)[1]]
  u <- tail(km$time[which(km$surv >= 1-p+e)],1)
  l <- km$time[which(km$surv <= 1-p-e)][1]
  S.u <- tail(km$surv[which(km$surv >= 1-p+e)],1)
  S.l <- km$surv[which(km$surv <= 1-p-e)][1]

  f <- (S.u-S.l) / (l-u)
  se <- se.S / f
  list("se"=se,"percentile"=perc)
}

# center
# treat ******
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

#Main: ######################################################################

  # Read in data:
  cgd <- read.csv("../../Data/cgd.csv")[,-c(1:2)]
  cat.i <- c(1,2,6,9)

  # Part 1: Logrank Test
  # KM Curve:
  km <- survfit(Surv(time,status) ~ treat, type="kaplan-meier",data=cgd)

  plot.km <- function(){
    survplot(km,col=c("purple","red"),lwd=1.2,lty=1,
         xlab="Time to Infection (Days)",ylab="Survival")
         title(main="KM Survival Curves for Infection Time")}

  # Find Median
  plac.i <- which(cgd$treat!="placebo")
  km.trmt <- survfit(Surv(time,status)~1, type="kaplan-meier",data=cgd[-plac.i,])
  km.plac <- survfit(Surv(time,status)~1, type="kaplan-meier",data=cgd[plac.i,])
  se.me.trmt <- se.pcnt(.5,km.trmt)
  se.me.plac <- se.pcnt(.5,km.plac) #NA

  CI.median.trmt <- qnorm(c(.025,.975),se.me.trmt$perc,se.me.trmt$se)
  km.trmt.median <- matrix(c(se.me.trmt$perc,CI.median.trmt),1,3)
  colnames(km.trmt.median) <- c("Estimate","CI.Lower","CI.Upper")

  # Logrank Test:
  logrank.treat <- survdiff(Surv(time,status) ~ treat, data=cgd)
  p.logrank.treat <- pchisq(q=logrank.treat$chisq,df=1,lower.tail=F)
  # p < .05 => treatment vs. placebo significantly different effects on survival.

  
  # Part 2: Cox Model
  # Model Selection: Stepwise Forward / Backward
  cox.ful <- coxph(Surv(time,status) ~ ., data=cgd)
  cox.red <- coxph(Surv(time,status) ~ 1, data=cgd)
  cox.stp <- step(cox.red,scope=list(lower=cox.red,upper=cox.ful),data=cgd,
                  direction="both")
  # Age is not significant. So, throw out.

  cox.mod <- coxph(Surv(time,status) ~ treat, data=cgd) # logrank test included
                                                        # same p value
  cox.mod.coef <- summary(cox.mod)$conf.int

  #predicted survival curve
  
  # Log.h ~ Log.T: Parallel => Proportional Hazards Assumption Met
  plot.log.H <- function(km,add=F,main="",col=1,type='p') {
    log.H <- log(-log(km$surv))
    log.T <- log(km$time)
    if (!add) {
      plot(log.T,log.H,col=col,pch=20,main=main,type=type,ylim=c(-5,1))
    } else {
      lines(log.T,log.H,col=col,pch=20,main=main,type=type)
    }
  }

  plot.model.ass <- function(){
    plot.log.H(km.trmt,col="red",type='s',
               main="Log Cumulative Hazard vs. Log Survival Time")
    plot.log.H(km.plac,col='purple',type='s',add=T)
    legend("topleft",col=c("red","purple"),lwd=3,legend=c("rIFN-g","Placebo"))
  }  

  #survplot(km,col=c("purple","red"),loglog=T,logt=T,add=T)
