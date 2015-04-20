# Read: P.203-214, 245-250
# Goal: Predict the balance of cardholders BEFORE issuing card
#       Determining characteristics of a cardholder that lead to high balances
#       Want: Moderate Balance, Avoid: High Balance, Don't Care about: Low Balance

rm(list=ls())
options("width"=120)
library(leaps)  #regsubsets

credit <- read.table("http://mheaton.byu.edu/Courses/Stat536/Case%20Studies/Credit/Data/Credit.csv",header=T,sep=",")
credit <- credit[,-1]
credit <- as.data.frame(credit)
#pairs(credit)

# A function to swap words separated by a colon
# the 'bigger' word goes before the colon
swap <- function(s){
  #p <- pos(":",s)
  p <- regexpr(":",s)
  if (p > 0 ){
    a   <- substr(s,1,p-1)
    b   <- substr(s,p+1,nchar(s))
    ifelse(a<b, paste(a,b,sep=":"), paste(b,a,sep=":")) 
  } else {
    s
  }
}

my.predict <- function(cof, dat){
  x <- model.matrix(Balance ~ .^2, data=dat)
  b <- cof
  colnames(x) <- sapply(colnames(x),swap); colnames(x) <- unname(colnames(x))
  names(b) <- sapply(names(b),swap)

  x <- x[, names(b)]

  pred <- x %*% b
  pred[pred < 0] <- 0
  
  pred
}

engine <- function(k){

  print(k)
  trainI <- sample(1:400,300)
  testI  <- setdiff(1:400,trainI)

  train <- credit[trainI,]
  test  <- credit[testI,]

  #Validate:
  # Using regsubsets():
  nvMax <- 25
  #crazy <- regsubsets(Balance ~ .^2, data=train, nvmax=nvMax, method="seqrep")    #3135,22 
  #crazy <- regsubsets(Balance ~ .^2, data=train, nvmax=nvMax, method="forward")   #3170,10 
  crazy <- regsubsets(Balance ~ .^2, data=train, nvmax=nvMax, method="backward")  #3244,15

  MSE <- NULL
  for (i in 1:nvMax){
    # my.predict function:
    cof <- coef(crazy,i)
    pred <- my.predict(cof,test)
    MSE[i] <- mean((test[,"Balance"]-pred)^2)
  }

  #par(mfrow=c(2,1))
  #plot(MSE,type='o',pch=20,xlab='Paramaters',main='MSE vs. Parameters')
  #plot((summary(crazy))$adjr2,type='o',pch=20,ylab=expression(paste(R^2)),xlab='Parameters',main=expression(paste('Adjusted ',R^2,' vs. Parameters')))
  min.MSE <- which.min(MSE)
  #min(which(abs(MSE-min(MSE)) < sd(MSE))) # the smallest number of parameters within 
  #                                        # one sd of the min MSE
  coefs <- coef(crazy,min.MSE)
  vars  <- names(coefs)

  vars
}

library(foreach)
library(doMC)
registerDoMC(16)

N <- 1000
vars <- foreach(n=1:N) %dopar% engine(n)

tab <- table(unlist(vars))
h <- N/2

my.plot <- function(){
  plot(tab,ylab="Count",las=2,cex.axis=.32,main="Parameter Counts Vs. Parameters")
  abline(h=h)
}

numParam <- sum(tab > h)
params <- names(which(tab > h))

terms <- paste(params[params!="(Intercept)"],collapse="+")
terms <- paste("Balance ~",terms)

trainI <- sample(1:400,300)
testI  <- setdiff(1:400,trainI)

train <- credit[trainI,]
test  <- credit[testI,]

X <- model.matrix(Balance ~ .^2, data=train)
X <- as.data.frame(cbind(train$Balance,X))
names(X)[1] <- "Balance"

Y <- model.matrix(Balance ~ .^2, data=test)
Y <- as.data.frame(cbind(test$Balance,Y))
names(Y)[1] <- "Balance"

mod <- lm(as.formula(terms),data=X)

pred <- predict(mod,Y,interval='prediction',level=.95)
mean((pred[,2] < Y$Balance) & (Y$Balance < pred[,3]))
