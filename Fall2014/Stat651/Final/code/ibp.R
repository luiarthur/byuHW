#rm(list=ls())
#library(expm)
library(MASS)
source("rfunctions.R")

tr <- function(m) sum(diag(m))

dmatnorm <- function(X,M,U,V){
  n <- nrow(X)
  p <- ncol(X)

  exp(-.5 * tr(solve(V) %*% t(X-M) %*% solve(U) %*% (X-M))) / 
  ( (2*pi)^(n*p) * det(V)^n * det(U)^p )^(.5)
}

ldmatnorm <- function(X,M,U,V){
  n <- nrow(X)
  p <- ncol(X)

  (-.5 * tr(solve(V) %*% t(X-M) %*% solve(U) %*% (X-M))) - 
  .5 * log((2*pi)^(n*p) * det(V)^n * det(U)^p)
}

#rmatnorm1 <- function(M,U,V){ # This Works! Not as fast.
#  n <- nrow(U)
#  p <- nrow(V)
#  Z <- matrix(rnorm(n*p),n,p)
#  X <- sqrtm(U) %*% Z %*% sqrtm(V) + M 
#  X
#}

rmatnorm <- function(M,U,V){ # Works. This is FASTER!
  # If M is nxp, then
  #    U is nxn
  #    V is pxp.

  vec.m <- cbind(as.vector(M))
  #v.x <- mvrnorm(1,vec(M),U %x% V)
  v.x <- vec.m+ t(chol(U %x% V))%*%rnorm(length(vec.m)) #MVN draw
  X <- matrix(v.x,nrow(U),ncol(V))
  X
}

#M <- matrix(1:20,4,5)
#U <- diag(4); V <- diag(5)
#N <- 10000
#X <- foreach(i=1:N) %dopar% rmatnorm(M,U,V) 
#(EX <- Reduce('+',X) / N)

r.X.given.Z <- function(Z,d=2,su=2,sv=2){ # X is n by d
  n <- nrow(Z)
  k <- ncol(Z)
  J <- matrix(1,k,d)
  Iu <- diag(n)
  Iv <- diag(d)

  rmatnorm(Z%*%J,su^2*Iu,sv^2*Iv)
} 

r.X.given.Z.given.a <- function(n=2,a=10,d=2,su=2,sv=2){ # X is n by d
  Z <- rIBP(n,a,F)
  k <- ncol(Z)
  #cat("\n"); print(paste("k = ", k))
  J <- matrix(1,k,d)
  Iu <- diag(n)
  Iv <- diag(d)

  rmatnorm(Z%*%J,su^2*Iu,sv^2*Iv)
}

count <- function(m,Zs){ # counts the appearance of the matrix m in a list
                         # of matrices Zs
  cnt <- 0
  for (i in 1:length(Zs)){
    if (all(dim(Zs[[i]]) == dim(m))){
      if (all(Zs[[i]] == m)) {
        cnt <- cnt + 1
      }
    }
  }
  cnt
}


normalize <- function(x) {
  sum <- 0
  for (i in 1:length(x)){
    sum <- sum + x[i]
  }  
  x / sum
}

find.k1 <- function(z){ # The number of new dishes sampled by each customer
  z <- as.matrix(z)
  k1 <- NULL
  k1[1] <- sum(z[1,])
  n <- nrow(z)
  max.col <- k1[1]
  for (i in 2:n){
     k1[i] <- sum(z[i,-c(1:max.col)]) # The number of new dishes sampled by cust. i
     if (k1[i]>0) max.col <- max.col+k1[i]
  }
  k1
}

find.kh <- function(z){

  # After finding kh,
  # foreach unique vector,
  # compute the factorial of the counts,
  # then multiply the quantities together.
  # So the constant in 14 is the product of the factorial of the counts of unique

  #n <- nrow(z)
  z <- as.matrix(z)
  k <- ncol(z)
  kh <- NULL
  z.col   <- list()
  for (i in 1:k) z.col[[i]] <- z[,i]
  uniq.zk <- unique(z.col)
  for (i in 1:length(uniq.zk)) kh[i] <- count(uniq.zk[[i]],z.col) 
  kh
}

dIBP <- function(z,a=1,log=F,exchangeable=T,const=T){
  z <- as.matrix(z[,which(apply(z,2,sum)>0)])
  n <- nrow(z)
  k <- ncol(z) # This is K+
  Hn <- sum(1/(1:n))

  if (k==0) {
    if (!log) {
      return(exp(-a*Hn))
    } else {
      return(-a*Hn)
    }  
  }

  mk <- apply(z,2,sum)
  kk <- NULL

  if (exchangeable) {
    kk <- find.kh(z) # exchangeable IBP uses kh
  } else {
    kk <- find.k1(z) # nonexchangeable uses k1
  }  

  A <- 0
  if (!log) {
    A <- ifelse(const,prod(gamma(kk+1)),1)
    density <- a^k / A * exp(-a*Hn) * prod(gamma(mk)*gamma(n-mk+1)/gamma(n+1))
  } else {
    A <- ifelse(const,sum(lgamma(kk+1)),0)
    density <- k*log(a) -A -a*Hn + sum(lgamma(mk) + lgamma(n-mk+1) - lgamma(n+1))
  }

  density
}


rIBP <- function(ppl=50,a=10,plotting=F,...){

  last <- rpois(1,a) # I assume that the first customer may sample NO dishes

  M <- matrix(0,ppl,last)
  M[1,0:last] <- 1 #V1

  for (n in 2:ppl){
    for (k in 0:last){
      mk <- sum(M[,k])
      M[n,k] <- sample(0:1,1,prob=c(1-mk/n,mk/n))
    }  
    newLast <- last+rpois(1,a/n)
    col0 <- matrix(0,ppl,newLast-last)
    if (ncol(col0) > 0){ # added & newLast > 0
      M <- cbind(M,col0)
      M[n,(last+1):newLast] <- 1
      last <- newLast
    }
  }

  if (plotting) {
    a.image(M,...)
    #pheatmap(M,cluster_row=F,cluster_col=F,display_numbers=T,number_format="%.0f")
  }

  M
}


#Z <- matrix(c(1,0,0,1,1,0,0,0,1),3,3)
#dIBP(Z,const=T)

#M <- matrix(1:12,3,4)
#U <- diag(3); V <- diag(4)
#N <- 10000
#
#X <- foreach(i=1:N) %dopar% rmatnorm(M,U,V)
#(EX <- Reduce('+',X) / N)

# Some assignment:
#  r.X.given.Z(n=5,a=12,d=3,su=1,sv=1); cat("\n\n")
#  plotGrid(10,1,plotting=F)

############################## Next Task: ###################################
#
# Changing Alpha:
#   1) E[X|Z=z]
#   2) E[X] = E[E[X|Z]] # all cells tend to a
#
# Metropolis:
#   1) X|Z ~ N(ZJ, s2_uI, s2I), where X(nxd), Z(nxk), J(kxd), I(dxd), I(kxk)
#   2) Z ~ IBP(a)
#
#   Q) What does X look like for various a, s2?
#   Q) What is Z|X?
#
#############################################################################


