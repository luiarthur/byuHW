sum.matrices <- function(Ms,return.matrices=F) { 
# Ms is a list of matrices of different lengths
# return.matrices is a boolean. If FALSE, function returns the sum of the matrices.
# If TRUE, function returns a list of the matrices also.

  l <- length(Ms)
  max.c <- max(unlist(lapply(Ms,ncol)))
  max.r <- max(unlist(lapply(Ms,nrow)))
  
  for (i in 1:l) {
    M <- Ms[[i]]
    
    ncol0 <- max.c - ncol(M)
    nrow0 <- max.r - nrow(M)

    if (ncol0>0) {
      col0 <- matrix(0,nrow(M),ncol0)
      M <- Ms[[i]] <- cbind(Ms[[i]],col0)
    }

    if (nrow0>0) {
      row0 <- matrix(0,nrow0,ncol(M))
      M <- Ms[[i]] <- rbind(Ms[[i]],row0)
    }
  }

  out <- Reduce("+",Ms)

  if (return.matrices) out <- list("sum"=out,"matrices"=Ms)

  out
}

art.roc <- function(y,p,n=1e3) {
  p0 <- seq(0,1,len=n)
  O <- matrix(0,n,2)
  colnames(O) <- c("1-Specificity","Sensitivity")

  for (i in 1:n) {
    p.hat <- ifelse(p>p0[i],1,0)

    true.ind <- which(y==1)
    false.ind <- which(y==0)
    sens <- mean(y[true.ind]==p.hat[true.ind])
    spec <- mean(y[false.ind]==p.hat[false.ind])

    O[i,1] <- 1-spec
    O[i,2] <- sens
  }
  O
}

library(pROC)
library(foreign)
dat_slc <- read.csv("data/slcmormon.txt")
dat_sf <- data.frame(read.spss("data/sfsample.SAV"))

### data cleaning
# initialize
x = data.frame(matrix(0, nrow(dat_slc) + nrow(dat_sf), 12))
colnames(x) = c("CITY", "LDS", "PRNTACTV", "SACRMTG", "AGE",
    "INCOME", "EDUC", "PRVPRAYR", "READSCRP", "MARITAL",
    "SEX", "FRIEND")
# 256, 28, 251, 252, 271, 36, 33, 270, 269, 3, 4, 6

# city indicator (0 for SF, 1 for SLC)
x$CITY = c(rep(0, nrow(dat_sf)), rep(1, nrow(dat_slc)))

# born in the church, or convert (7)
temp_sf = ifelse(dat_sf[,7] == "Born in the LDS Church, but no longer affiliated" |
    dat_sf[,7] == "Life-long member and still affiliated", "birth", "convert")
temp_slc = ifelse(dat_slc[,7] == 1 | dat_slc[,7] == 2, "birth", "convert")
x$LDS = c(temp_sf, temp_slc)

# parents activity level (256)
x$PRNTACTV = c(dat_sf[,256], dat_slc[,256])

# sacrament meeting attendance (28)
x$SACRMTG = c(dat_sf[,28], dat_slc[,28])

# age, continu-ized (251)
x$AGE = c(dat_sf[,251], dat_slc[,251])
for (i in 1:length(x$AGE)){
    if (is.na(x$AGE[i]))
        next
    if (x$AGE[i] == 1) x$AGE[i] = 20
    if (x$AGE[i] == 2) x$AGE[i] = 25
    if (x$AGE[i] == 3) x$AGE[i] = 30
    if (x$AGE[i] == 4) x$AGE[i] = 35
    if (x$AGE[i] == 5) x$AGE[i] = 40
    if (x$AGE[i] == 6) x$AGE[i] = 45
    if (x$AGE[i] == 7) x$AGE[i] = 50
    if (x$AGE[i] == 8) x$AGE[i] = 55
    if (x$AGE[i] == 9) x$AGE[i] = 60
    if (x$AGE[i] == 10) x$AGE[i] = 65
    }

# income continu-ized (252) (smaller)
x$INCOME = c(dat_sf[,252], dat_slc[,252])
for (i in 1:length(x$INCOME)){
    if (is.na(x$INCOME[i]))
        next
    if (x$INCOME[i] == 1) x$INCOME[i] = 4000
    if (x$INCOME[i] == 2) x$INCOME[i] = 5000
    if (x$INCOME[i] == 3) x$INCOME[i] = 6000
    if (x$INCOME[i] == 4) x$INCOME[i] = 7000
    if (x$INCOME[i] == 5) x$INCOME[i] = 8000
    if (x$INCOME[i] == 6) x$INCOME[i] = 9000
    if (x$INCOME[i] == 7) x$INCOME[i] = 10000
    if (x$INCOME[i] == 8) x$INCOME[i] = 11000
    if (x$INCOME[i] == 9) x$INCOME[i] = 15000
    if (x$INCOME[i] == 10) x$INCOME[i] = 20000
    if (x$INCOME[i] == 11) x$INCOME[i] = 35000
    if (x$INCOME[i] == 12) x$INCOME[i] = 50000
    }

# education level (271)
x$EDUC = c(dat_sf[,271], dat_slc[,271])

# private prayer (36)
x$PRVPRAYR = c(dat_sf[,36], dat_slc[,36])

# scripture reading (33)
x$READSCRP = c(dat_sf[,33], dat_slc[,33])

# marital status (270)
x$MARITAL = c(dat_sf[,270], dat_slc[,270])

# sex (269)
x$SEX = 2 - c(dat_sf[,269], dat_slc[,269])

# priesthood office (3)
# x$PSTH = c(dat_sf[,3], dat_slc[,3])
# x$PSTH[x$SEX == 0] = 0  # set priesthood to 0 if female
#                         # may need to deal with this
#                         # just ignore priesthood office?

# lds friends (4)
x$FRIEND = c(c(dat_sf[,4])-1, dat_slc[,4])

# lds friends (6)
# x$MISSION = c(dat_sf[,6], dat_slc[,6]) - 1

# response
temp_sf = ifelse(dat_sf[,7] == "Born in the LDS Church, but no longer affiliated" |
    dat_sf[,7] == "Once converted, but no longer affiliated", 0, 1)
temp_slc = ifelse(dat_slc[,7] == 1 | dat_slc[,7] == 4, 0, 1)
y = c(temp_sf, temp_slc)

# mising values
miss.y = which(is.na(y))
miss.x = NULL
for (i in 1:ncol(x))
    miss.x = unique(c(miss.x, which(is.na(x[,i]))))
# every row that has at least one missing value
miss.all = c(miss.y, miss.x)

# counts for responses after removing missing values
table(y[-miss.all])

# check to make sure they are removed
sum(is.na(y[-miss.all]))
sum(is.na(x[-miss.all,]))

y = y[-miss.all]
x = x[-miss.all,]

# make all the x's factors (except age and income)
x[,names(x)] = lapply(x[,names(x)], factor)
x$AGE = c(x$AGE)
x$INCOME = c(x$INCOME)
x$FRIEND = c(x$FRIEND)


### model fitting
#mod = glm(y ~ ., data = x, family = binomial(link=logit))
#summary(mod)
#
#mod = glm(y ~ CITY + SEX + AGE + INCOME, data = x, family = binomial)
#summary(mod)

fit.model <- function(my.link="logit") {
  mod.int <- glm(y~1, data=x, family=binomial(link=my.link))
  mod.all <- glm(y~., data=x, family=binomial(link=my.link))
  mod.forward <- step(mod.int, scope=list(lower=mod.int,upper=mod.all), 
                      direction='both',data=x) 
  #summary(mod)

  resp <- predict(mod.forward,type="response")
  my.roc <- roc(y,resp)
  #plot(my.roc)
  my.auc <- auc(my.roc)+0
  
  out <- list("model"=mod.forward,"roc"=my.roc,"auc"=round(my.auc,5))
}


mod.logit <- fit.model("logit") 
mod.probit <- fit.model("probit") 
mod.cloglog <- fit.model("cloglog") 

summary(mod.logit$mod)

par(mfrow=c(3,1))
  plot(mod.logit$roc,main=paste("Logit Model,AUC =",mod.logit$auc),col="blue")
  legend("right",legend=names(mod.logit$model$coef),cex=.9,title="Covariates")
  plot(mod.probit$roc,main=paste("Probit Model,AUC =",mod.probit$auc),col="blue")
  legend("right",legend=names(mod.probit$model$coef),cex=.8,title="Covariates")
  plot(mod.cloglog$roc,main=paste("Clog-log Model,AUC =",mod.cloglog$auc),col="blue")
  legend("right",legend=names(mod.cloglog$model$coef),cex=.75,title="Covariates")
par(mfrow=c(1,1))

# Notes from Arthur:
#   - Should we go with the Logit model because AUC is just slightly (<.01) smaller
#     than clog-log and probit models, but is easier to interpret, and has the 
#     fewest number of covariates?
##############################################################################



final.mod <- function(my.link="logit") {
  temp <- glm(y~CITY+PRVPRAYR+CITY*PRVPRAYR+FRIEND+SACRMTG+LDS,data=x,
              family=binomial(link=my.link))
  resp <- predict(temp,type="response")
  my.roc <- roc(y,resp)
  #plot(my.roc)
  my.auc <- auc(my.roc)+0

  out <- list("model"=temp,"roc"=my.roc,"auc"=round(my.auc,5))
}

my.mod.logit <- final.mod()
my.mod.probit <- final.mod("probit")
my.mod.cloglog <- final.mod("cloglog")

summary(my.mod.logit$mod)
summary(my.mod.probit$mod)
summary(my.mod.cloglog$mod)

library(xtable)
xtable(my.mod.logit$mod)

make.my.plots <- function() {
  p.log <- art.roc(y,predict(my.mod.logit$mod,type="response"))
  plot(p.log,main=paste("Logit Model, AUC =",my.mod.logit$auc),col="blue",
       type="l",lwd=3)
  abline(0,1)     
}
 
pdf("auc.pdf")
  make.my.plots() 
dev.off()

pdf("auc3.pdf",height=3.5)
  par(mfrow=c(1,3))

  p.log <- art.roc(y,predict(my.mod.logit$mod,type="response"))
  plot(p.log,main=paste("Logit Model, AUC =",my.mod.logit$auc),col="blue",
       type="l",lwd=3)
  abline(0,1)

  p.pro <- art.roc(y,predict(my.mod.probit$mod,type="response"))
  plot(p.pro,main=paste("Probit Model, AUC =",my.mod.probit$auc),col="blue",
       type="l",lwd=3)
  abline(0,1)

  p.clog <- art.roc(y,predict(my.mod.cloglog$mod,type="response"))
  plot(p.clog,main=paste("Cloglog Model, AUC =",my.mod.cloglog$auc),col="blue",
       type="l",lwd=3)
  abline(0,1)

  par(mfrow=c(1,1))
dev.off()


cv <- function(train.ind,dat=x) { # train.ind = index for training set
  train <- dat[train.ind,]
  test <- dat[-train.ind,]
  y.train <- y[train.ind]
  y.test <- y[-train.ind]
  
  mod <- glm(y.train~CITY+PRVPRAYR+CITY*PRVPRAYR+FRIEND+SACRMTG+LDS,data=train,
              family=binomial(link="logit"))
  
  pred <- predict(mod,newdata=test,type="response")

  my.roc <- roc(y.test,pred)
  my.auc <- auc(my.roc)+0
  artROC <- art.roc(y.test,pred)


  out <- list("model"=mod,"roc"=my.roc,"auc"=round(my.auc,5),"artROC"=artROC,
              "pred"=pred,"y.test"=y.test)
  out
}

n <- nrow(x)
####################################

samp <- sample(1:n)
l <- 10
k <- floor(n/l)
cv.mods <- as.list(1:l)
for (i in 1:l) {
  samp.i <- samp[((i-1)*k+1):(min(n,i*k))]
  print(samp.i)
  cv.mods[[i]] <- cv(samp.i)
}

errorRates <- NULL
library(doMC)
library(foreach)
registerDoMC(16)
do.it <- function(i,p=.5) {
  errorRates[i] <- mean(cv.mods[[i]]$y.test == ifelse(cv.mods[[i]]$pred>p,1,0))
}

errorRates <- foreach(i=1:l,.combine=rbind,.errorhandling="pass") %dopar% do.it(i,.5)
mean.error <- mean(errorRates)
mean.aucs <- mean(unlist(lapply(cv.mods,function(x) x$auc)))
sd.aucs <- sd(unlist(lapply(cv.mods,function(x) x$auc)))

range(errorRates)
range(aucs)

plot(cv.mods[[3]]$artROC,main=paste("Logit Model (cv), AUC =",cv.mods[[3]]$auc),
     col="blue", type="l",lwd=3); abline(0,1)
####################################### MICKEY
m = 500
lower = 0
upper = 1
cutoff = seq(lower,upper,length=m)
pos.rate = double(m)
neg.rate = double(m)
error = double(m)

p <- predict(my.mod.logit$mod,type="response")
for (i in 1:m){
    t.pos = sum(y == 1 & p > cutoff[i])
    t.neg = sum(y == 0 & p < cutoff[i])
    total.neg = (n-0)-sum(y)
    total.pos = sum(y)

    pos.rate[i] = 1-t.pos/total.pos
    neg.rate[i] = 1-t.neg/total.neg
    error[i] = 1-(t.pos+t.neg)/(n-0)
    }

pdf("error.pdf")
plot(cutoff, pos.rate, type='s', col="darkblue", ylab="Error Rate",
    xlab="Threshold", main="Predictive Error Rates")
points(cutoff, neg.rate, type='s', col="red")
points(cutoff, error, type='s', col="green", lwd=2)
legend(0.6, 1, legend=c("False Positive", "False Negative", "Overall Error"),
    col=c("darkblue", "red", "green"), lwd=c(1,1,2), lty=1, cex=1.1)
dev.off()


#####################################################3
cutoff = .5
p <- predict(my.mod.logit$mod,type="response")
t.pos = sum(y == 1 & p > cutoff)
t.neg = sum(y == 0 & p < cutoff)
total.neg = (n-0)-sum(y)
total.pos = sum(y)

pos.rate[i] = 1-t.pos/total.pos
neg.rate[i] = 1-t.neg/total.neg
error[i] = 1-(t.pos+t.neg)/(n-0)

t.pos
t.neg

confusion <- matrix(0,2,2)
colnames(confusion) <- c("Prediction=0","Prediction=1")
rownames(confusion) <- c("y=0","y=1")

confusion[1,1] <- t.neg
confusion[1,2] <- total.neg-t.neg
confusion[2,1] <- total.pos-t.pos
confusion[2,2] <- t.pos

confusion
xtable(confusion)
