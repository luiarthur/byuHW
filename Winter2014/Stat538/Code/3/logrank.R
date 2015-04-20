library(survival)

#### FAKE DATA 

g1 <- c(3.1,6.8,9,9,11.3,16.2)
c1 <- c(1,0,1,1,0,1)

g2 <- c(8.7,9,10.1,12.1,18.7,23.1)
c2 <- c(1,1,0,0,1,0)

x <- c(rep(1,6),rep(0,6))
time <- c(g1,g2)
status <- c(c1,c2)

### informal comparison: KM curves

g1.KM <- survfit(Surv(g1,c1)~1,type="kaplan-meier")
g2.KM <- survfit(Surv(g2,c2)~1,type="kaplan-meier")
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

#postscript(file="aml_KM_compare.ps",height=10,width=10)
plot(g1.KM.time,g1.KM.surv,ylim=c(ymin,ymax),type="s",lty=1,
xlab="Time",
ylab="Survival",
main="Example",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)
lines(g1.KM.time,g1.KM.lower,type="s",lty=2)
lines(g1.KM.time,g1.KM.upper,type="s",lty=2)
lines(g2.KM.time,g2.KM.surv,lty=1,type="s",col=2)
lines(g2.KM.time,g2.KM.lower,type="s",lty=2,col=2)
lines(g2.KM.time,g2.KM.upper,type="s",lty=2,col=2)
#dev.off()


### formal comparison of survival curves using log-rank test
res.logrank <- survdiff(Surv(time,status)~x)
res.logrank

### formal comparison of survival curves using Peto-Peto generalized Wilcoxon test
res.peto <- survdiff(Surv(time,status)~x,rho=1)
res.peto


#### UIS DATA

# treat=0 is the short program and treat=1 is the long program
# The patients were randomly assigned to two different sites (site=0 is site A and site=1 is site B).
# age indicates age at enrollment,
# herco indicates heroine or cocaine use in the past three months (herco=1 indicates heroine and cocaine use, herco=2 indicates either heroine or cocaine use and herco=3 indicates neither heroine nor cocaine use)
# ndrugtx indicates the number of previous drug treatments.

uis.dat <- read.csv("uis.csv")
head(uis.dat)

g1.KM <- survfit(Surv(uis.dat$time[uis.dat$treat==0],uis.dat$censor[uis.dat$treat==0])~1,type="kaplan-meier")
g2.KM <- survfit(Surv(uis.dat$time[uis.dat$treat==1],uis.dat$censor[uis.dat$treat==1])~1,type="kaplan-meier")
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

max(g1.KM.time)
max(g2.KM.time)

#postscript(file="aml_KM_compare.ps",height=10,width=10)
plot(g2.KM.time,g2.KM.surv,ylim=c(ymin,ymax),type="s",lty=1,
xlab="Time",
ylab="Survival",
main="Example",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)
#lines(g2.KM.time,g2.KM.lower,type="s",lty=2)
#lines(g2.KM.time,g2.KM.upper,type="s",lty=2)
lines(g1.KM.time,g1.KM.surv,lty=1,type="s",col="red")
#lines(g1.KM.time,g1.KM.lower,type="s",lty=2,col=2)
#lines(g1.KM.time,g1.KM.upper,type="s",lty=2,col=2)
##dev.off()
### red = treat0 = short program
### black = treat1 = long program

### formal comparison of survival curves using log-rank test
res.logrank <- survdiff(Surv(uis.dat$time,uis.dat$censor)~uis.dat$treat)
res.logrank

### formal comparison of survival curves using Peto-Peto test
res.logrank <- survdiff(Surv(uis.dat$time,uis.dat$censor)~uis.dat$treat, rho=1)
res.logrank

