nurse <- read.csv("nursing.csv")

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
#nurse.full <- coxph(Surv(duration,completion) ~ .^2, data=nurse,x=T)
#nurse.1    <- coxph(Surv(duration,completion) ~   1, data=nurse,x=T)
#step.mod   <- step(nurse.1,scope=list(lower=nurse.1,upper=nurse.full),data=nurse,direction="both")


nurse.mod <- coxph(Surv(duration,completion) ~ race+poverty+smoked+
                                               education+poverty*education,
                                               data=nurse,x=T)
# poverty:education not significant
nurse.red <- coxph(Surv(duration,completion) ~ race+poverty+smoked+education,data=nurse,x=T)

test.stat <- -2*(nurse.red$loglik[2] - nurse.full$loglik[2])
p.val <- pchisq(test.stat,df=length(nurse.full$coefficients)-
                length(nurse.red$coefficients),lower.tail=F)


# PLOTS:
# Poverty:  
plot.pov <- function(cex){
  pov.comp <- survfit(Surv(duration,completion) ~ poverty,
                      type="kaplan-meier",data=nurse)
  plot(pov.comp,col=2:3,lwd=2,main="Survival for Mothers in Porverty",
       xlab="Time (Weeks)",ylab="Survival")
  legend("topright",legend=c("Mother Not in Poverty","Mother in Poverty"),
         col=2:3,lwd=3,cex=cex)
}

# Race:  
plot.race <- function(cex){  
  race.comp <- survfit(Surv(duration,completion) ~ race,
                       type="kaplan-meier",data=nurse)
  plot(race.comp,col=2:4,lwd=2,main="Survival for Different Races",
       xlab="Time (Weeks)",ylab="Survival")
  legend("topright",legend=c("White","Black","Other"),cex=cex,col=2:4,lwd=3)
}

# Smoked:  
plot.smoked <- function(cex){  
  smoked.comp <- survfit(Surv(duration,completion) ~ smoked,
                         type="kaplan-meier",data=nurse)
  plot(smoked.comp,col=2:3,lwd=2,main="Survival for Smoking Mothers",
       xlab="Time (Weeks)",ylab="Survival")
  legend("topright",legend=c("Mother Not Smoking",
                             "Mother Smoking"),
                             cex=cex,col=2:4,lwd=3)
}

# Education:
plot.edu <- function(cex){  
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

