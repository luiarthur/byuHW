rmvn <- function(n=1,mu=0,Sigma=1){
  draws <- mu + t(chol(Sigma)) %*% rnorm(n)
  draws
}

dingam <- function(x,a,b)
  b^a/gamma(a) * x^(-a-1) * exp(-b/x)

ringam <- function(n,a,b)
  1/rgamma(n,a,b)

n <- 1000
b0 <- 20
b1 <- 3
s2 <- 7
beta <- c(b0,b1)
e <- rnorm(n,sd=sqrt(s2))
X <- cbind(1,rnorm(n,.01))
Y <- X%*%beta + e
plot(X[,2],Y,col='red',pch=20)

mod <- lm(Y~X[,-1])
abline(mod,col="blue")

m <- c(0,0) # Prior mean for beta
S <- matrix(c(1,0,0,1),2,2)
a <- 2
b <- 2 # beta is a scale parameter, not rate!
n <- length(Y)

XTX <- t(X)%*%X
XTY <- t(X)%*%Y

update.beta <- function(s2,m,S){
  Si <- solve(S)
  m.new <- solve(XTX / s2 + Si) %*% (XTY / s2 + Si%*%m)
  S.new <- solve(XTX / s2 + Si)
  beta.new <- rmvn(2,m.new,S.new)
  list("m"=m.new,"S"=S.new,"b"=beta.new)
}

update.s2 <- function(B,a,b) {
  a.new <- n/2 + a
  b.new <- t(Y-X%*%B) %*% (Y-X%*%B)/2 + b
  s2.new <- ringam(1,a.new,b.new) #1/rgamma(1,a,b) 
  list("a"=a.new,"b"=b.new,"s2"=s2.new)
}

gibbs <- function(a,b,m,S,B,s2,N=50000){
  out <- matrix(0,N,3)
  for (i in 2:N){

    up.beta <- update.beta(s2,m,S)
    m <- up.beta$m
    S <- up.beta$S
    beta <- up.beta$b

    up.s2   <- update.s2(beta,a,b)
    a <- up.s2$a
    b <- up.s2$b
    s2 <- up.s2$s2

    out[i,] <- c(beta[1],beta[2],s2)
  }
  list("B"=beta,"m"=m,"S"=S,"s2"=s2,"a"=a,"b"=b,"sampler"=out)
}

out <- gibbs(a,b,m,S,c(0,0),1)
abline(out$B[1],out$B[2],col='gold',lwd=3)
result <- rbind(c(out$B,out$s2),c(beta,s2),c(mod$coef,summary(mod)$sigma^2))
colnames(result) <- c("b0","b1","s2")
rownames(result) <- c("Bayes","Truth","MLE")
result
