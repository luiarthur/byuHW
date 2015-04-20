set.seed(5)

n <- 15
k <- 3
N <- n*k
K <- 2*k
sig.r2 <- .1

gam <- c(-1,-3,4,0,2,-2)
#gam <- c(gam,-sum(gam))
b <- c(1,1)
G <- diag(36,K)

#X <- cbind(1,matrix(rnorm(N,5,.5),N))
X <- cbind(1,matrix(c(runif(n,0,1),runif(n,1,2),runif(n,2,3)),N))

Zz <- matrix(0,N,k)
for (k in 1:k) {
  Zz[((k-1)*(N/3)+1):(k*(N/3)),k] <- 1
}

Zx <- matrix(0,N,k)
for (kk in 1:k) {
  Zx[((kk-1)*(N/3)+1):(kk*(N/3)),kk] <- X[((kk-1)*(N/3)+1):(kk*(N/3)),2]
}

Z <- cbind(Zz,Zx)

e <- rnorm(n,0,sqrt(sig.r2))

y <- X%*%b + Z%*%gam + e

R <- diag(1,N)
V <- Z%*%G%*%t(Z) + R

b.hat <- solve(t(X) %*% solve(V) %*% X) %*% t(X) %*% solve(V) %*% y
gam.hat <- G %*% t(Z) %*% solve(V) %*% (y-X%*%b.hat)

y.hat <- X%*%b.hat + Z%*%gam.hat 
plot(y-y.hat)

#abline(lm(y~X[,2]),col="red")
clust.num <- apply(Zz,1,which.max) 
plot(X[,2],y,col=clust.num,pch=20)
#K <- ncol(Z)
for (kk in 1:k) {
  ind <- which(clust.num==kk)
  abline(b.hat[1]+Z[ind,1:k]%*%gam.hat[1:k], b.hat[2]+gam.hat[k+kk],col=kk,lwd=2)
}

cbind(b,b.hat)
cbind(gam,gam.hat)


