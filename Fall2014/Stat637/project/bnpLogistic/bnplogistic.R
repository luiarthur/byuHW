source("functions/plotpost.R")
source("data/gendata.R")

#library(ElemStatLearn)
#data(ozone)
#y <- ifelse(ozone$o>70,1,0)
#X <- ozone[,-1]
#pairs(ozone)
#pairs(cbind(y,X))

dat <- gendat(100)
b <- dat$b
y <- dat$y
x <- dat$x
X <- cbind(1,x)
plot(x,y)


