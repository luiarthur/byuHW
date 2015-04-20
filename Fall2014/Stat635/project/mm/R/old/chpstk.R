source("rfunctions.R")
source("gibbs2.R")

dat <- read.table("chopstick.txt")
y <- dat[,1]
X <- cbind(1,dat[,2])

B <- 500
elapsed.time <- system.time(out<-gibbs.post(y,X,B=B,showProgress=T,plotProgress=T))

EZ <- est.Z(out$Zs)
EZ <- clust.Z(EZ)
a.image(EZ,axis.num=F,main="Posterior Estimate of Z")

G <- diag(100,ncol(EZ))
R <- diag(mean(out$sig.r2),length(y))
V <- EZ %*% G %*% t(EZ) + R
beta.hat <- solve(t(X) %*%solve(V) %*%X)%*%t(X) %*%solve(V) %*%y
gam.hat <- G%*%t(EZ)%*%solve(V)%*%(y-X%*%beta.hat)

clust.num <- apply(EZ,1,which.max) 
plot(X[,2],y)
plot(X[,2],y,col=clust.num)
K <- ncol(EZ)
for (kk in 1:K) {
  ind <- which(clust.num==kk)
  abline(beta.hat[1]+EZ[ind,]%*%gam.hat, beta.hat[2],col=kk)
}

beta.hat
gam.hat
