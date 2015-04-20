library(xtable)

y1 <- sleep$extra[which(sleep$group==1)]
y2 <- sleep$extra[which(sleep$group==2)]

y1 <- ifelse(y1>0,1,0)
y2 <- ifelse(y2>0,1,0)

n <- length(y1)

#1a) 
calc.ci <- function(p,sd) {
  qnorm(c(.025,.975),p,sd)
}

# p1:
    p1.hat <- mean(y1)
    y1 <- sum(y1)
  #i)
    sd1..5 <- sqrt(.5*(1-.5) / n)
  #ii)
    sd1.phat <- sqrt(p1.hat*(1-p1.hat) / n)
  #iii)
    ptil1 <- (y1+2)/(n+4) 
    sd1.ptil <- sqrt(ptil1*(1-ptil1)/(n+4))
  #iv)
    eta1.hat <- log(y1/(n-y1))
    eta1 <- qnorm(c(.025,.957),eta1.hat,sqrt(1/y1 + 1/(n-y1)))
    p1.eta.ci <- exp(eta1) / (1+exp(eta1))
  #Compare  
    X1 <- matrix(0,4,4)
    rownames(X1) <- c("p=.5","p.hat","p.til","eta")
    colnames(X1) <- c("Estimate","SD","95%.CI.Lower","95%.CI.Upper")
    X1[,1] <- c(.5,p1.hat,ptil1,eta1.hat)
    X1[,2] <- c(sd1..5,sd1.phat,sd1.ptil,NA)
    X1[1:3,3:4] <- t(apply(X1[1:3,],1,function(x) calc.ci(x[1],x[2])))
    X1[4,3:4] <- p1.eta.ci
    X1

# p2:
    p2.hat <- mean(y2)
    y2 <- sum(y2)
  #i)
    sd2..5 <- sqrt(.5*(1-.5) / n)
  #ii)
    sd2.phat <- sqrt(p2.hat*(1-p2.hat) / n)
  #iii)
    ptil2 <- (y2+2)/(n+4) 
    sd2.ptil <- sqrt(ptil2*(1-ptil2)/(n+4))
  #iv)
    eta2.hat <- log(y2/(n-y2))
    eta2 <- qnorm(c(.025,.957),eta2.hat,sqrt(1/y2 + 1/(n-y2)))
    p2.eta.ci <- exp(eta2) / (1+exp(eta2))
  #Compare  
    X2 <- matrix(0,4,4)
    rownames(X2) <- c("p=.5","p.hat","p.til","eta")
    colnames(X2) <- c("Estimate","SD","95%.CI.Lower","95%.CI.Upper")
    X2[,1] <- c(.5,p2.hat,ptil2,eta2.hat)
    X2[,2] <- c(sd2..5,sd2.phat,sd2.ptil,NA)
    X2[1:3,3:4] <- t(apply(X2[1:3,],1,function(x) calc.ci(x[1],x[2])))
    X2[4,3:4] <- p2.eta.ci
    X2
  
#1b) Comment on how these ci compare.
#1c) Do any of the ci types provide wvidence that one of the drugs is significantly
#    better?

#3)
dat <- read.table("case.txt",header=T)
y <- dat$Rem
m <- dat$Cases
x <- dat$L1

#a)
mod.clog <- glm(y/m~x,weights=m,family=binomial(link="cloglog"))
mod.log  <- glm(y/m~x,weights=m,family=binomial(link="logit"))
mod.prob <- glm(y/m~x,weights=m,family=binomial(link="probit"))

# logit has slightly higher resid dev. than probit. But logit is easier to 
# interpret, so I choose logit.

#b)
smod.log <- summary(mod.log)
smod.log
# Interpret: The expected increase in log odds is .14486 for a one level increase
#            in L1. The coefficient is significant for predicting probability of 
#            remission.
new.dat <- data.frame(x=seq(26,26.1,len=1000))
pred <- predict(mod.log,newdat=new.dat,type="response")
l1.5 <- new.dat[which.min(abs(pred-.5)),]
l1.5


pred.8 <- predict(mod.log,newdat=data.frame(x=8),se.fit=T)
est <- pred.8$fit
est.se <- pred.8$se.fit
eta.ci <- calc.ci(est,est.se)
pred.8.ci <- exp(eta.ci) / (1+exp(eta.ci))
pred.8.ci
