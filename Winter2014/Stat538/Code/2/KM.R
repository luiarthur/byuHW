library(survival)

################################
# no censoring
################################

time.example <- c(2,2,3,5,5,7,9,16,16,18)
cens.example <- c(1,1,1,1,1,1,1,1,1,1)

# creation of survival object
Surv(time.example,cens.example)

# creation of KM curve based on survival object
survfit(Surv(time.example,cens.example)~1,type="kaplan-meier")

summary(survfit(Surv(time.example,cens.example)~1,type="kaplan-meier"))

KM.example <- summary(survfit(Surv(time.example,cens.example)~1,type="kaplan-meier"))

ymax <- 1
ymin <- 0

KM.example.time <- c(0,KM.example$time)
KM.example.surv <- c(1,KM.example$surv)

postscript(file="ch4_KM_example_nocens.ps",height=10,width=10)
plot(KM.example.time,KM.example.surv,ylim=c(ymin,ymax),type="s",lty=1,
xlab="Time",
ylab="Survival",
main="K-M Survival Curve (no censoring)",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)
dev.off()


################################
# censoring
################################

time.example <- c(2,2,3,5,5,7,9,16,16,18)
cens.example <- c(1,1,0,1,0,1,1,1,1,0)

Surv(time.example,cens.example)
survfit(Surv(time.example,cens.example)~1,type="kaplan-meier")
summary(survfit(Surv(time.example,cens.example)~1,type="kaplan-meier"))

KM.example <- summary(survfit(Surv(time.example,cens.example)~1,type="kaplan-meier"))

ymax <- 1
ymin <- 0

KM.example.time <- c(0,KM.example$time)
KM.example.surv <- c(1,KM.example$surv)

postscript(file="ch4_KM_example_cens.ps",height=10,width=10)
plot(KM.example.time,KM.example.surv,ylim=c(ymin,ymax),type="s",lty=1,
xlab="Time",
ylab="Survival",
main="K-M Survival Curve",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)
dev.off()



################################
# censoring: UIS data set
# The goal of the UIS  data, and of the smaller version called uis_small that we are using here, is to model time
#      until return to drug use for patients enrolled in two different residential treatment programs that
#      differed in length (treat=0 is the short program and treat=1 is the long program).  The patients were
#      randomly assigned to two different sites (site=0 is site A and site=1 is site B).  The variable age
#      indicates age at enrollment, herco indicates heroine or cocaine use in the past three months (herco=1
#      indicates heroine and cocaine use, herco=2 indicates either heroine or cocaine use and herco=3 indicates
#      neither heroine nor cocaine use) and ndrugtx indicates the number of previous drug treatments.  The
#      variable time contains the time until return to drug use and the censor variable indicates whether the
#      subject returned to drug use (censor=1 indicates return to drug use and censor=0 otherwise).
################################

uis.dat <- read.csv("uis.csv")
head(uis.dat)
uis.relapse <- uis.dat$time
uis.cens <- uis.dat$censor

summary(survfit(Surv(uis.relapse,uis.cens)~1,type="kaplan-meier"))
KM.example <- summary(survfit(Surv(uis.relapse,uis.cens)~1,type="kaplan-meier"))

ymax <- 1
ymin <- 0

KM.example.time <- c(0,KM.example$time)
KM.example.surv <- c(1,KM.example$surv)

postscript(file="KM.uis.ps",height=10,width=10)
plot(KM.example.time,KM.example.surv,ylim=c(ymin,ymax),type="s",lty=1,
xlab="Time",
ylab="Survival",
main="K-M Survival Curve",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)
dev.off()


################################
# censoring with CIs
################################

uis.dat <- read.csv("uis.csv")
head(uis.dat)
uis.relapse <- uis.dat$time
uis.cens <- uis.dat$censor

summary(survfit(Surv(uis.relapse,uis.cens)~1,type="kaplan-meier"))
KM.example <- summary(survfit(Surv(uis.relapse,uis.cens)~1,type="kaplan-meier"))

ymax <- 1
ymin <- 0

KM.example.time <- c(0,KM.example$time)
KM.example.surv <- c(1,KM.example$surv)

KM.example.upper <- c(1,KM.example$upper)
KM.example.lower <- c(1,KM.example$lower)

postscript(file="KM.uis.ps",height=10,width=10)
plot(KM.example.time,KM.example.surv,ylim=c(ymin,ymax),type="s",lty=1,
xlab="Time",
ylab="Survival",
main="K-M Survival Curve",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)
lines(KM.example.time,KM.example.lower,type="s",lty=2)
lines(KM.example.time,KM.example.upper,type="s",lty=2)
dev.off()


################################
# median survival time
################################

Surv(uis.relapse,uis.cens)
survfit(Surv(uis.relapse,uis.cens)~1,type="kaplan-meier") # median=166
summary(survfit(Surv(uis.relapse,uis.cens)~1,type="kaplan-meier")) # median=166 (smallest surv time for which the surv function is less than or equal to 0.50)

# largest surv time for which surv function is greater than or equal to 0.50+0.01=0.51 = 161
# smallest surv time for which surv function is less than or equal to 0.50-0.01=0.49 = 168

# f.hat = (S(161) - S(168))/(168-161) = (0.51 - 0.489)/7 = 0.021/7 = 0.003
# se(mp) = (1/f.hat) * se(S(mp) = (0.01995)/.003 = 6.65

# 95% CI: 166 +/- 1.96*6.65 = 166 +/- 13.034 = c(153,179)
# narrower than 95% CI reported by R (c(148,184))

################################
# hazard function
################################

uis.dat <- read.csv("uis.csv")
head(uis.dat)
uis.relapse <- uis.dat$time
uis.cens <- uis.dat$censor

KM.uis <- summary(survfit(Surv(uis.relapse,uis.cens)~1,type="kaplan-meier",conf.type="log-log"))

time.minus1 <- c(0,KM.uis$time[-length(KM.uis$time)])
time.hazard <- KM.uis$time
interval.width <- time.hazard - time.minus1
inst.hazard <- KM.uis$n.event/((KM.uis$n.risk)*(interval.width))
length(KM.uis$surv)
length(inst.hazard)

postscript("ch4_KM_example_hazard.ps",height=10,width=10)
plot(time.hazard,inst.hazard,type="s",
xlab="t",
ylab="h(t)",
main="Hazard Function",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)
dev.off()


library(muhaz)
hazard.uis <- kphaz.fit(uis.relapse,uis.cens)
kphaz.plot(hazard.uis)


################################
# hazard function (constant)
################################

time.example <- c(1,1,2,3,4,5,6,7,8,9)
cens.example <- c(1,1,0,0,0,1,0,0,0,0)

Surv(time.example,cens.example)
survfit(Surv(time.example,cens.example)~1,type="kaplan-meier")
summary(survfit(Surv(time.example,cens.example)~1,type="kaplan-meier"))

KM.example <- summary(survfit(Surv(time.example,cens.example)~1,type="kaplan-meier",conf.type="log-log"))

time.minus1 <- c(0,KM.example$time)
time.hazard <- KM.example$time
interval.width <- time.hazard - time.minus1[-length(time.minus1)]
inst.hazard <- KM.example$n.event/((KM.example$n.risk)*(interval.width))

postscript("ch4_KM_example_hazard.ps",height=10,width=10)
plot(time.hazard,inst.hazard,type="s",
xlab="t",
ylab="h(t)",
main="Hazard Function",
cex.main=2,
cex.axis=1.5,
cex.lab=1.5)
dev.off()

