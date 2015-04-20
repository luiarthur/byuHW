########################################################################
# attempt to fit Cox model using large data set

library(survival)
dlbcl.data <- as.matrix(read.table("dlbcl_expressiondata.csv",sep=","))
dlbcl.surv <- as.matrix(read.table("dlbcl_survdata.csv",sep=","))
dim(dlbcl.data)
dlbcl.surv

coxph(Surv(dlbcl.surv[,1],dlbcl.surv[,2])~1)
coxph(Surv(dlbcl.surv[,1],dlbcl.surv[,2])~t(dlbcl.data[1:5,]))
coxph(Surv(dlbcl.surv[,1],dlbcl.surv[,2])~t(dlbcl.data[1:150,]))
coxph(Surv(dlbcl.surv[,1],dlbcl.surv[,2])~t(dlbcl.data[1:300,]))
coxph(Surv(dlbcl.surv[,1],dlbcl.surv[,2])~t(dlbcl.data[1:300,]),coxph.control(iter.max=100))


###

# generate data from two identical populations

n <- 1000
reject.res <- p.value <- NA
for (i in 1:10000) {
  x <- rnorm(n,5,3)
  y <- rnorm(n,5,3)
  p.value[i] <- t.test(x,y,"two.sided")$p.value
  reject.res[i] <- ifelse(p.value[i] < .05,1,0)
}
sum(reject.res)/length(reject.res) # % declared as significant


###

# experiment-wide significance level
# calculate familywise error rate

n.tests <- 10
alpha <- .05
1-(1-alpha)^n.tests # probability of Type I error

n.tests <- 100
alpha <- .05
1-(1-alpha)^n.tests # probability of Type I error

n.tests <- 1000
alpha <- .05
1-(1-alpha)^n.tests # probability of Type I error

###

# Bonferroni correction

alpha <- .05
n.tests <- 10
alpha/n.tests

alpha <- .05
n.tests <- 1000
alpha/n.tests

alpha <- .05
n.tests <- 100000
alpha/n.tests

