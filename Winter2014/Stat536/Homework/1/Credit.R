# Read: P.203-214, 245-250
# Goal: Predict the balance of cardholders BEFORE issuing card
#       Determining characteristics of a cardholder that lead to high balances
#       Want: Moderate Balance, Avoid: High Balance, Don't Care about: Low Balance

rm(list=ls())
options("width"=120)
library(car)    #vif
library(leaps)  #regsubsets

credit <- read.table("http://mheaton.byu.edu/Courses/Stat536/Case%20Studies/Credit/Data/Credit.csv",header=T,sep=",")
credit <- credit[,-1]
credit <- as.data.frame(credit)
#credit$Cards <- paste(credit$Cards)
#credit$Cards[credit$Cards>="7"] <- "7+"
#head(credit)
#pairs(credit)
#
#Response:
#Balance (Continuous)
#
#Input:
#  Quantitative: |  Qualitative: 
#    Income      |   
#    Limit       |    Gender
#    Rating      |    Student
#    Age         |    Married
#    Education   |    Ethnicity
#table(credit$Cards)
#table(credit$Gender)
#table(credit$Student)
#table(credit$Married)
#table(credit$Ethnicity)

trainI <- sample(1:400,300)
testI  <- setdiff(1:400,trainI)

train <- credit[trainI,]
test  <- credit[testI,]

train$Gender <- factor(train$Gender)
train$Student <- factor(train$Student)
train$Married <- factor(train$Married)
train$Ethnicity <- factor(train$Ethnicity)

par(mfrow=c(2,1))
lower.mod <- lm(Balance ~ 1, data=train)
#upper.mod <- lm(Balance ~ ., data=train)
#
#summary(upper.mod)
#
#forwardS <- step(lower.mod,scope=list(lower=lower.mod,upper=upper.mod),upper=upper.mod,direction="both")
#mod <- eval(forwardS$call)
#
#summary(mod)
#
#vif(mod)
#  # Limit  15.44
#  # Rating 15.47
#
#hist(resid(mod)) #Residuals not symmetric, looks Gamma
#qqnorm(resid(mod))
# How do I get the Bias?
# How do I get the MSE?
 
# Throw out Rating because of high VIC?

# Interaction Model
int.mod <- lm(Balance ~ .^2, data=train) 
#forwardS <- step(lower.mod,scope=list(lower=lower.mod,upper=int.mod),upper=int.mod,direction="both",k=log(400)) #BIC
forwardS <- step(lower.mod,scope=list(lower=lower.mod,upper=int.mod),upper=int.mod,direction="both",k=log(400)) #BIC
#forwardS <- step(lower.mod,scope=list(lower=lower.mod,upper=int.mod),upper=int.mod,direction="both") #AIC

mod <- eval(forwardS$call)
summary(mod)
hist(resid(mod))
qqnorm(resid(mod))
par(mfrow=c(1,1))

# Test the test set for mod2
#pos <- function(c,s){
#  p <- 0
#  while( (substr(s,p,p) != c) & (p < nchar(s) )){
#    p <- p + 1
#  }
#  ifelse(substr(s,p,p)==c,p,0)
#}

#Validate:
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

x <- model.matrix(Balance ~ .^2, data=test)
b <- mod$coef
colnames(x) <- sapply(colnames(x),swap); colnames(x) <- unname(colnames(x))
names(b) <- sapply(names(b),swap)

x <- x[, names(b)]

pred <- x %*% b
pred[pred<0] <- 0
MSE <- mean((test[,"Balance"]-pred)^2)

# Using regsubsets():
nvMax <- 25
crazy <- regsubsets(Balance ~ .^2, data=train, nvmax=nvMax, method="seqrep")
MSE <- NULL

for (i in 1:nvMax){
  cof <- coef(crazy,i)

  x <- model.matrix(Balance ~ .^2, data=test)
  b <- cof 
  colnames(x) <- sapply(colnames(x),swap); colnames(x) <- unname(colnames(x))
  names(b) <- sapply(names(b),swap)

  x <- x[, names(b)]

  pred <- x %*% b
  pred[pred<0] <- 0

  MSE[i] <- mean((test[,"Balance"]-pred)^2)
}

MSE

#terms <- paste(names(cof)[-1],collapse="+")
#terms <- paste("Balance ~",terms)
#lm(as.formula(terms),data=train)
