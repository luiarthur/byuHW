source("plotPost.R")
source("countDown.R")
dat <- read.table("sore.txt",header=T)

y <- dat$Y
X <- cbind(1,dat$D,dat$T)

m.probit <- function(y,X,n=nrow(X),k=ncol(X),prior.b=cbind(rep(0,k),rep(1,k)),
                      cand.s=rep(.01,k),B=1e4,burn=round(B*.1),
                      trim.burn=F) {
  
  calc.dev <- function(p) -2*sum(y*log(p)+(1-y)*log(1-p))

  out <- matrix(0,B,k)
  log.g <- function(b,m,s,p) {
    sum(y*log(p)+(1-y)*log(1-p)) - ((b-m)/s)^2/2
  }
  acc <- rep(0,k)
  dev <- NULL
  dev[1] <- calc.dev(.5)
  p <- NULL

  for (i in 2:B) {
    old.time <- Sys.time()
    out[i,] <- out[i-1,]
    for (j in 1:k) {
      cand <- rnorm(1,out[i,j],cand.s[j])
      cand.b <- out[i,]
      cand.b[j] <- cand
      p.cand <- pnorm(c(X%*%cand.b))
      p <- p.prev <- pnorm(c(X%*%out[i,]))
      m.ratio <- log.g(cand,prior.b[j,1],prior.b[j,2],p.cand) - 
                 log.g(out[i,j],prior.b[j,1],prior.b[j,2],p.prev)
      
      if (!(any(p.cand==0) || any(p.cand==1))) {
        if (m.ratio > log(runif(1))) {
          out[i,j] <- cand
          p <- p.cand
          if (!(trim.burn)) {
            acc[j] <- acc[j] + 1
          } else if (i>burn){
            acc[j] <- acc[j] + 1
          }
        }
      }

    }

    dev[i] <- calc.dev(p)
    count.down(old.time,i,B)
  }

  acc.rate <- acc/B
  if (trim.burn) acc.rate <- acc/(B-burn)
  list("post"=out,"acc"=acc.rate,"dev"=dev)
}

#X1 <- X[,2] - mean(X[,2])
#c.X <- cbind(X[,1],X1,X[,3])
#out <- m.probit(y,c.X,cand.s=c(1,.1,1))
out <- m.probit(y,X,cand.s=c(1,.02,1),trim.burn=T,B=1e4)
out$acc

#opts <- par(no.readonly=T)
#plot.posts(out$post,names=c("b0","b1","b2"));par(opts)

#hpd.95 <- t(apply(out$post,2,get.hpd))
#rownames(hpd.95) <- paste0("beta",0:2)
#colnames(hpd.95) <- c("Lower 95% HPD","Upper 95% HPD")
#hpd.95

# Need to compute deviance!!!
# Need to interpret values!!!

dev <- out$dev[which(!(is.na(out$dev)))]
#plot.post(dev);par(opts)

x0 <- c(1,44,1)
p <- pnorm(out$post%*%x0)
#plot.post(p,"p");par(opts)

pdf("pics/postMs.pdf")
  plot.posts(out$post,names=c("b0","b1","b2"),cex.l=.6,cex.a=.8)
dev.off()
pdf("pics/postMDev.pdf")
  plot.post(dev)
dev.off()
pdf("pics/postMP.pdf")
  plot.post(p,"p")
dev.off()
