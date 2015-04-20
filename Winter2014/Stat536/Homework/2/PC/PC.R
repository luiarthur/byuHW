rm(list=ls())
options("width"=150)
library(pls)

GDP <- read.csv("../GDP_data.csv",sep=",")
gdp <- GDP[,-c(1:2)]

X <- model.matrix(GR6096 ~ .,data=gdp)[,-1]
x.bin.ind <- which(apply(X,2,function(x) length(unique(x))) <= 2)
X.bin <- X[, x.bin.ind] # Binary X
X.qnt <- X[,-x.bin.ind] # Quantitative X
Y <- gdp$GR6096

# R Functions:
#pcr.mod <- pcr(Y~X.qnt,scale=F,validation="CV")
#validationplot(pcr.mod,val.type="MSEP")
#
#pcr.mod.2 <- pcr(Y~X.qnt,scale=T,ncomp=11)
#validationplot(pcr.mod.2,val.type="MSEP")
#summary(pcr.mod.2)

n <- nrow(X)
P <- ncol(X)
M <- 11
psi <- t(eigen(cov(X))$vectors[,1:M])
Sig <- 0
mean.X <- apply(X,2,mean)

for (i in 1:n){
  Sig <- Sig + (X[i,]-mean.X) %*% t(X[i,]-mean.X)
}
Sig <- Sig/(n-1)
X <- scale(X)

Z <- matrix(0,n,M)
for (i in 1:n){
  for (m in 1:M){
    for (p in 1:P){
      Z[i,m] <- Z[i,m] + psi[m,p]*X[i,p]   
    }
  }
}

Y <- scale(Y)
mod <- lm(Y~Z)
theta <- coef(mod)

beta <- rep(0,P+1)
beta[1] <- theta[1]
for (p in 2:length(beta)){
  for (m in 2:length(theta)){
    beta[p] <- theta[m] * psi[m-1,p-1]
  }
}

beta <- matrix(beta,ncol=1)
rownames(beta) <- c("(Intercept)",colnames(X))
head(cbind(beta[rev(order(abs(beta))),]),11)
# beta is now the reg. coef. for the scaled and centered X's and Y's.
