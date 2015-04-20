set.seed(5)
gendata <- function(n.vec=c(30,30,30),g=c(-7,1,6),b=c(5,4),
                    G=diag(36,3),x=rnorm(sum(n.vec),5,.5),sigr=1,
                    e=rnorm(sum(n.vec),0,sigr),R=diag(sigr^2,sum(n.vec))) {

  n <- sum(n.vec)
  Z <- matrix(0,n,3)
  nv <- cumsum(c(0,n.vec))
  for (k in 1:length(g)) {
    Z[(nv[k]+1):nv[k+1],k] <- 1
  }
  
  X <- cbind(1,x)
  y <- X%*%b + Z%*%g + e
  V <- Z%*%G%*%t(Z) + R

  list("y"=y,"X"=X,"b"=b,"Z"=Z,"g"=g,"sigr2"=sigr^2)
}

plot.mm <- function(y,x,b,z,g,addplot=F,line=T,...) {
  clust.num <- apply(z,1,which.max) 

  if (addplot) {
    points(x,y,col=clust.num+1,...)
  } else {
    plot(x,y,col=clust.num+1,...)
  }

  if (line) {
    K <- ncol(z)
    for (kk in 1:K) {
      ind <- which(clust.num==kk)
      abline(b[1]+z[ind,]%*%g, b[2],col=kk+1,lwd=2)
    }
    abline(lm(y~x),lwd=2,col="grey80")
  }
}


