#Using the data and identified model(s) from Homework4,
#assess model fit and appropriateness through residual analysis.
#Analysis should include:
#1) assessment of the proportional hazards assumption, (cox.zph)
#2) assessment of overallfit,
#3) identification and discussion of influential observations.

nurse <- read.csv("../4/nursing.csv")

library(survival)
#pairs(nurse)

nurse <- nurse[,-8]
nurse$race <- as.factor(nurse$race)
nurse$poverty <- as.factor(nurse$poverty)
nurse$smoked <- as.factor(nurse$smoked)
nurse$alcohol <- as.factor(nurse$alcohol)
nurse$prenatal <- as.factor(nurse$prenatal)

#nurse.full <- coxph(Surv(duration,completion) ~ .,data=nurse,x=T)
# poverty:education kinda significant

# Stepwise: Forward / Backward selection:
nurse.full <- coxph(Surv(duration,completion) ~ .^2, data=nurse,x=T)
nurse.1    <- coxph(Surv(duration,completion) ~   1, data=nurse,x=T)
#step.mod   <- step(nurse.1,scope=list(lower=nurse.1,upper=nurse.full),data=nurse,direction="both")

nurse.mod <- coxph(Surv(duration,completion) ~ race+poverty+smoked+
                                               education+poverty*education,
                                               data=nurse,x=T)
# poverty:education not significant 
# THIS WILL BECOME OUR MODEL:
nurse.red <- coxph(Surv(duration,completion) ~ race+poverty+smoked+education,data=nurse,x=T)

#step.stat <- -2*(nurse.red$loglik[2] - step.mod$loglik[2])
#step.p.val <- pchisq(temp.stat,df=length(step.mod$coefficients)-
#                length(nurse.red$coefficients),lower.tail=F)

test.stat <- -2*(nurse.red$loglik[2] - nurse.full$loglik[2])
p.val <- pchisq(test.stat,df=length(nurse.full$coefficients)-
                length(nurse.red$coefficients),lower.tail=F)


# PLOTS:
# Poverty:  
plot.pov <- function(cex=1){
  pov.comp <- survfit(Surv(duration,completion) ~ poverty,
                      type="kaplan-meier",data=nurse)
  plot(pov.comp,col=2:3,lwd=2,main="Survival for Mothers in Porverty",
       xlab="Time (Weeks)",ylab="Survival")
  legend("topright",legend=c("Mother Not in Poverty","Mother in Poverty"),
         col=2:3,lwd=3,cex=cex)
}

# Race:  
plot.race <- function(cex=1){  
  race.comp <- survfit(Surv(duration,completion) ~ race,
                       type="kaplan-meier",data=nurse)
  plot(race.comp,col=2:4,lwd=2,main="Survival for Different Races",
       xlab="Time (Weeks)",ylab="Survival")
  legend("topright",legend=c("White","Black","Other"),cex=cex,col=2:4,lwd=3)
}

# Smoked:  
plot.smoked <- function(cex=1){  
  smoked.comp <- survfit(Surv(duration,completion) ~ smoked,
                         type="kaplan-meier",data=nurse)
  plot(smoked.comp,col=2:3,lwd=2,main="Survival for Smoking Mothers",
       xlab="Time (Weeks)",ylab="Survival")
  legend("topright",legend=c("Mother Not Smoking",
                             "Mother Smoking"),
                             cex=cex,col=2:4,lwd=3)
}

# Education:
plot.edu <- function(cex=1){  
  edu <- as.factor((nurse$education > 12) * 1)
  #edu <- nurse$education
  education.comp <- survfit(Surv(duration,completion) ~ edu,
                            type="kaplan-meier",data=nurse)
  plot(education.comp,col=2:3,lwd=2,main="Survival for Mothers (Education)",
       xlab="Time (Weeks)",ylab="Survival")
  legend("topright",legend=c("Mothers <= 12 Years of Education",
                             "Mothers >  12 Years of Education"),
                              col=2:4,lwd=3,cex=cex)
}




# PART2:

  # RESIDUALS:####################################################################

  # Estimation of Baseline Hazard and Baseline Survival:
  # summary(survfit(nurse.red)) # baseline hazard
  # coxph.detail(nurse.red)$hazard # hazard increment

  
  #1) Test the Proportional Hazards Assumption for a Cox Regression Model Fit
  edu <- as.factor((nurse$education > 12) * 1)
  cox.zph(nurse.red)
  # Since education is significant, dichotomize, then stratify
  nurse.red.dic <- coxph(Surv(duration,completion) ~ race+poverty+smoked+edu,
                                                     data=nurse,x=T)
  cox.zph(nurse.red.dic)


  nurse.red.dic.strat <- coxph(Surv(duration,completion) ~ 
                               race+poverty+smoked+strata(edu),data=nurse,x=T)
  cox.zph(nurse.red.dic.strat)
  #################################################################################

  #2) Assessment of Overallfit: using deviance residuals
  ########## Martingale residuals
  #nurse.red$residuals
  #residuals(nurse.red,type="martingale") # default
  plot.mart <- function(){
    plot(nurse.red$residuals,ylab="Martingale Residuals",pch=20)
  }

  ########## Deviance residuals
  #residuals(nurse.red,type="deviance")
  plot.dev <- function(){
    plot(residuals(nurse.red,type="deviance"),ylab="Deviance Residuals",pch=20)
    identify(residuals(nurse.red,type="deviance"))
  }

  

  #3) identification and discussion of influential observations: 
      # using identify
      # using scores
  getlikdis <- function(x) {
    likdis <- NA
    for (i in 1:nrow(residuals(x,type="score"))) {
      likdis[i] <- t(as.vector(residuals(x,type="score")[i,])) %*% x$var %*%
        as.vector(residuals(x,type="score")[i,])
    }
    return(likdis)
  }

  #as.vector(residuals(nurse.red,type="score")[1,])
  #t(as.vector(residuals(nurse.red,type="score")[1,])) %*% nurse.red$var %*%
  #  as.vector(residuals(nurse.red,type="score")[1,])

  likdis.nurse.red <- getlikdis(nurse.red) #need this
  plot(likdis.nurse.red,xlab="Time",ylab="Likelihood Displacement",
       main="Likelihood Displacement Vs. Time") # need this
  abline(h=.15)
  nurse.new <- as.data.frame(nurse[likdis.nurse.red < .15,])

  edu.new <- as.factor((nurse.new$education > 12) * 1)
  nurse.red.new <- coxph(Surv(duration,completion) ~ race+poverty+smoked+
                                                 strata(edu.new),data=nurse.new,x=T)
  nurse.red.new
  nurse.red.dic.strat

  likdis.nurse.new <- getlikdis(nurse.red.new)
  plot(likdis.nurse.new,xlab="Time",ylab="Likelihood Displacement")
  abline(h=.15)

