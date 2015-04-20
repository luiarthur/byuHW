library(truncnorm) # for rtruncnorm
source("plotPost.R")
source("countDown.R")

dat <- read.table("sore.txt",header=T)

y <- dat$Y
X <- cbind(1,dat$D,dat$T)

mvrnorm <- function(M,S,n=nrow(S)) M + t(chol(S)) %*% rnorm(n)

updateZ <- function(x,y,b){
  z <- numeric(length(y))
  z[y==1] <- rtruncnorm(sum(y==1),a=0,b=Inf,mean=x[y==1,] %*% b,sd=1)
  z[y==0] <- rtruncnorm(sum(y==0),a=-Inf,b=0,mean=x[y==0,] %*% b,sd=1)
  z
}


gibb <- function(y,X,n=nrow(X),k=ncol(X),B=1e4,burn=round(B*.1),
                 trim.burn=F,V=diag(1,k)) {
  

  calc.dev <- function(p) -2*sum(y*log(p)+(1-y)*log(1-p))
  S <- solve(t(X)%*%X + solve(V))
  Xt <- t(X)

  # Initialize Parameters
  z <- 1
  beta <- matrix(0,B,k)
  dev <- NULL
  dev[1] <- calc.dev(.5)
  p <- NULL
  #######################

  for (i in 2:B){
    #Updates:
    old.time <- Sys.time()
    z <- updateZ(X,y,beta[i-1,])
    beta[i,] <- mvrnorm(S %*% Xt%*%z, S) 
    dev[i] <- calc.dev(pnorm(X%*%beta[i,]))
    count.down(old.time,i,B)
  }

  list("beta"=beta,"dev"=dev)
}

outs <- gibb(y,X,B=1e4)
out <- outs$beta

#opts <- par(no.readonly=T)
#plot.posts(out,names=c("b0","b1","b2"),cex.a=1,tck.d=2); par(opts)

#hpd.95 <- t(apply(out,2,get.hpd))
#rownames(hpd.95) <- paste0("beta",0:2)
#colnames(hpd.95) <- c("Lower 95% HPD","Upper 95% HPD")
#hpd.95

dev <- outs$dev[which(!(is.na(outs$dev)))]
#plot.post(dev); par(opts)

x0 <- c(1,44,1)
z <- out%*%x0 
p <- pnorm(z)

#plot.post(p,"p");par(opts)

source("plotPost.R")
pdf("pics/postGibbs.pdf")
  plot.posts(out,names=c("b0","b1","b2"),cex.a=.8,tck.d=2,cex.l=.6)
dev.off()
pdf("pics/postGibbsDev.pdf")
  plot.post(dev)
dev.off()
pdf("pics/postGibbsP.pdf")
  plot.post(p,"p")
dev.off()

