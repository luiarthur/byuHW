library(lme4)
library(MASS)
set.seed(5)

n <- 3*30
# THIS SET WORKS with sig.r2=.1!!!
#X <- cbind(1,matrix(runif(n,0,5),n))
#b <- c(2,4)
#gam <- c(3,7,30)

#b <- c(2,4)
gam <- c(-7,1,6)
b <- c(5,4)

#gam.vec <- cbind(gam[1]+rnorm(1000),gam[2]+rnorm(1000),gam[3]+rnorm(1000))
#
#sum.x <- matrix(0,3,3)
#mean.x <- apply(gam.vec,2,mean)
#for (i in 1:(n/3)) {
#  sum.x <- sum.x + (gam.vec[i,]) %*% t(gam.vec[i,])
#}
#S <- 1/(n/3-1) * sum.x
#G <- S
G <- diag(36,3)
#apply(mvrnorm(1000,c(0,0,0),G),2,mean)

X <- cbind(1,matrix(rnorm(n,5,.5),n))
#X <- matrix(runif(n,0,5),n)

Z <- matrix(0,n,3)
for (k in 1:length(gam)) {
  Z[((k-1)*(n/3)+1):(k*(n/3)),k] <- 1
}

e <- rnorm(n,0,1)

y <- X%*%b + Z%*%gam + e

R <- diag(1,n)
V <- Z%*%G%*%t(Z) + R

b.hat <- solve(t(X) %*% solve(V) %*% X) %*% t(X) %*% solve(V) %*% y
gam.hat <- G %*% t(Z) %*% solve(V) %*% (y-X%*%b.hat)

y.hat <- X%*%b.hat + Z%*%gam.hat 
plot(y-y.hat)

plot(X[,2],y,pch=20,main="y=Xb+e",xlab="x")
abline(lm(y~X[,2]),lwd=2)

#abline(lm(y~X[,2]),col="red")
clust.num <- apply(Z,1,which.max) 
plot(X[,2],y,col=clust.num+1,pch=20,main="y=Xb+Zg+e",xlab="x")
K <- ncol(Z)
for (kk in 1:K) {
  ind <- which(clust.num==kk)
  abline(b.hat[1]+Z[ind,]%*%gam.hat, b.hat[2],col=kk+1,lwd=2)
}
abline(b.hat[1],b.hat[2],lwd=2)

cbind(b,b.hat)
cbind(gam,gam.hat)

#ind(y,X[,2],apply(Z,1,which.max))
#sink("data.txt")
#  cbind(y,X[,2],apply(Z,1,which.max))
#sink()
