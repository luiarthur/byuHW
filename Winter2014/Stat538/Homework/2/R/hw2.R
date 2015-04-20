#The final write up should be in either pdf or Word format and should include both a
#brief introduction (describing the data) and a summary, consisting of comments on
#the results.

rm(list=ls())
library(survival)
options("width"=120)
system("rm -f *.pdf")

mel <- read.csv("../melanoma.csv")

bcg <- mel[which(mel$Treatment.Received==1),]
cbp <- mel[which(mel$Treatment.Received==2),]

#par(mfrow=c(2,1))

# Resmission Duration:
   km.bcg.1 <- summary(survfit(Surv(bcg$Remission.Duration,bcg$rcens)~1,type="kaplan-meier"))
   km.cbp.1 <- summary(survfit(Surv(cbp$Remission.Duration,cbp$rcens)~1,type="kaplan-meier"))
  km1.bcg.med <- km.bcg.1$time[which(km.bcg.1$surv < .5)[1]]
  km1.cbp.med <- km.cbp.1$time[which(km.cbp.1$surv < .5)[1]]

  # Plotting:
  pdf("remission.pdf")
    xmax <- max(km.bcg.1$time, km.cbp.1$time)
    plot (c(0,km.bcg.1$time), c(1,km.bcg.1$surv), ylim=0:1, xlim=c(0,xmax),
          type="s", xlab="Remission Duration (months)", ylab="Survival", 
          main="K.M. Survival Curve for Remission Duration", col="blue",lwd=2)
    lines(c(0,km.bcg.1$time),c(1,km.bcg.1$lower),type="s",col="blue",lwd=2,lty=3)
    lines(c(0,km.bcg.1$time),c(1,km.bcg.1$upper),type="s",col="blue",lwd=2,lty=3)

    lines(c(0,km.cbp.1$time),c(1,km.cbp.1$surv),type="s",col="red",lwd=2)
    lines(c(0,km.cbp.1$time),c(1,km.cbp.1$lower),type="s",col="red",lwd=2,lty=3)
    lines(c(0,km.cbp.1$time),c(1,km.cbp.1$upper),type="s",col="red",lwd=2,lty=3)

    legend("bottomleft",legend=c("Bacillus Calmette-Guerin","Corynebacterium parvum",
           "Confidence Interval"), col=c("blue","red","black"),lwd=3,lty=c(1,1,3))
  dev.off()

# Survival Time:
  km.bcg.2 <- summary(survfit(Surv(bcg$Survival.Time,bcg$survcens)~1,type="kaplan-meier"))
  km.cbp.2 <- summary(survfit(Surv(cbp$Survival.Time,cbp$survcens)~1,type="kaplan-meier"))

  km2.bcg.med <- km.bcg.2$time[which(km.bcg.2$surv < .5)[1]]
  km2.cbp.med <- km.cbp.2$time[which(km.cbp.2$surv < .5)[1]]

  # Plotting:
  pdf("survival.pdf")
    xmax <- max(km.bcg.2$time, km.cbp.2$time)
    plot (c(0,km.bcg.2$time), c(1,km.bcg.2$surv), ylim=0:1, xlim=c(0,xmax),
          type="s", xlab="Survival Time (months)", ylab="Survival", 
          main="K.M. Survival Curve for Survival Time", col="blue",lwd=2)
    lines(c(0,km.bcg.2$time),c(1,km.bcg.2$lower),type="s",col="blue",lwd=2,lty=3)
    lines(c(0,km.bcg.2$time),c(1,km.bcg.2$upper),type="s",col="blue",lwd=2,lty=3)

    lines(c(0,km.cbp.2$time),c(1,km.cbp.2$surv),type="s",col="red",lwd=2)
    lines(c(0,km.cbp.2$time),c(1,km.cbp.2$lower),type="s",col="red",lwd=2,lty=3)
    lines(c(0,km.cbp.2$time),c(1,km.cbp.2$upper),type="s",col="red",lwd=2,lty=3)

    legend("bottomleft",legend=c("Bacillus Calmette-Guerin","Corynebacterium parvum",
           "Confidence Interval"), col=c("blue","red","black"),lwd=3,lty=c(1,1,3))
  dev.off()

  km.median <- c(km1.bcg.med,km1.cbp.med,km2.bcg.med,km2.cbp.med)
  # median: 6.4, 15.9, 19.5,  ___
  CI.cbp.1 <- 15.9 + c(-1,1) * 1.96 * 5.2785
  # CI:     ___, (5.55414, 26.24586), ____, ___ 


  se.pcnt <- function(p,km,e=.05){# p is a percentage for a percentile
    # Computes the standard error of a given percentile
    se.S <- km$std.err[which(km$surv<p)[1]]
    u <- tail(km$time[which(km$surv >= 1-p+e)],1)
    l <- km$time[which(km$surv <= 1-p-e)][1]
    S.u <- tail(km$surv[which(km$surv >= 1-p+e)],1)
    S.l <- km$surv[which(km$surv <= 1-p-e)][1]

    f <- (S.u-S.l) / (l-u)
    se <- se.S / f
    se
  } # This works!!! :)
  #se <- se.pcnt(.5,km.cbp.1)
  #CI.cbp.1 <- qnorm(c(.025,.975),km.median[2],se)
  

# Fisher Test:  
  mel$Initial.Stage <- substr(mel$Initial.Stage,1,1)
  mel$Age <- (mel$Age < mean(mel$Age)) * 1

  fisher.test(table(mel$Treatment.Received,mel$Initial.Stage))
  fisher.test(table(mel$Treatment.Received,mel$Age))
  fisher.test(table(mel$Treatment.Received,mel$Gender))

# Since the fisher test does not return significant p-values,
# there's no need to stratify out data.
# logrank:
  remission.logrank <- survdiff(Surv(mel$Remission.Duration,mel$rcens) ~ 
                                mel$Treatment.Received)
  p.remission.logrank <- pchisq(q=remission.logrank$chisq,df=1,lower.tail=F)

  survival.logrank <- survdiff(Surv(mel$Survival.Time,mel$survcens) ~
                               mel$Treatment.Received)
  p.survival.logrank <- pchisq(q=survival.logrank$chisq,df=1,lower.tail=F)

