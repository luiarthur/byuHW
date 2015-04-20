M <- matrix(c(132,176,127,
              42,6,12,
              172,129,130,
              56,4,15),4,3,byrow=T)

rownames(M) <- c("MW","MB","FW","FB")
colnames(M) <- c("D","R","I")

n <- sum(M)
y <- numeric(n)
X <- matrix(0,n,3)
X[,1] <- 1

rM <- nrow(M)
cM <- ncol(M)
k <- 0
for (i in 1:rM) {
  gender <- substr(rownames(M)[i],1,1)
  race <- substr(rownames(M)[i],2,2)
  for (j in 1:cM) {
    party <- colnames(M)[j]

    current.counter <- 0
    while (current.counter < M[i,j]) {
      current.counter <- current.counter + 1
      k <- k + 1
      y[k] <- party
      X[k,2] <- race
      X[k,3] <- gender
    }  
  }
}

dat <- cbind(y,X)

library(nnet)
mod <- multinom(y~X[,2]+X[,3])
summary(mod)

#Call:
#multinom(formula = y ~ X[, 2] + X[, 3])
#
#Coefficients:
#  (Intercept)  X[, 2]W   X[, 3]M
#I   -1.388246 1.118285 0.2201885
#R   -2.565342 2.278131 0.5727615
#
#Std. Errors:
#  (Intercept)   X[, 2]W   X[, 3]M
#I   0.2296432 0.2335159 0.1582525
#R   0.3436657 0.3427926 0.1575206
#
#Residual Deviance: 2085.782 
#AIC: 2097.782 


