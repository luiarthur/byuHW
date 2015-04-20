source("rfunctions.R")
source("ibp.R")
source("gibbs.R")

source("data/combine.R",chdir=T)
elapsed.time <- system.time(out <- gibbs.post(Y,a=1,B=5000,burn=0,showProgress=T,
                                              plotProgress=T,a.a=3,a.b=2,
                                              siga=1,sigx=.5))

# What Next?

M <- out$Zs
alpha <- out$alpha
burn <- 1000#round(length(M) * .1)

n.col <- unlist(lapply(M,ncol))

pdf("draw.post.out/trace.pdf")
par(mfrow=c(1,2))
  plot(n.col,type="l",main="Trace Plot: Number of Columns in Z",lwd=1,cex=.1,
       col="gray30",pch=20)
  mean.col <- round(mean(n.col[-(1:burn)]),4)
  var.col <-  round(var(n.col[-(1:burn)]),5)
  abline(h=mean.col,lwd=2,col="red")
  legend("bottomright",legend=c(paste("Mean=",mean.col), 
                                paste("Variance =",var.col)),title.col="gray30",
                                title=paste("After Burn-in of",burn,":"),bty="n")

  plot.post(alpha[-(1:burn)])
  #plot(alpha,type="l",main="Trace Plot: Alpha",lwd=1,cex=.1,
  #     col="gray30",pch=20)
  #mean.a <- round(mean(alpha[-(1:burn)]),4)
  #var.a <-  round(var(alpha[-(1:burn)]),5)
  #abline(h=mean.a,lwd=2,col="red")
  #legend("topleft",legend=c(paste("Mean=",mean.a),
  #                              paste("Variance =",var.a)),title.col="gray30",
  #                              title=paste("After Burn-in of",burn,":"),bty="n")
par(mfrow=c(1,1))
dev.off()

pdf("draw.post.out/traceplot.pdf")
  plot(n.col,type="l",main="Trace Plot: Number of Columns in Z",lwd=1,cex=.1,
       col="gray30",pch=20)
  mean.col <- round(mean(n.col[-(1:burn)]),4)
  var.col <-  round(var(n.col[-(1:burn)]),5)
  abline(h=mean.col,lwd=2,col="red")
  legend("bottomright",legend=c(paste("Mean=",mean.col), 
                                paste("Variance =",var.col)),title.col="gray30",
                                title=paste("After Burn-in of",burn,":"),bty="n")
dev.off()

pdf("draw.post.out/tracealpha.pdf")
  plot(alpha,type="l",main="Trace Plot: Alpha",lwd=1,cex=.1,
       col="gray30",pch=20)
  mean.a <- round(mean(alpha[-(1:burn)]),4)
  var.a <-  round(var(alpha[-(1:burn)]),5)
  abline(h=mean.a,lwd=2,col="red")
  legend("topleft",legend=c(paste("Mean=",mean.a),
                                paste("Variance =",var.a)),title.col="gray30",
                                title=paste("After Burn-in of",burn,":"),bty="n")
dev.off()

EAXZ <- function(X,Z,siga=1,sigx=.5) {
  k <- ncol(Z)
  Ik <- diag(k)
  ZT <- t(Z)
  out <- solve(ZT%*%Z +(sigx/siga)^2 * Ik,ZT%*%X)
  #out <- solve(ZT%*%Z, ZT%*%X)
  out
}

Z.post <- M[-(1:burn)] # Burn in about 100
Z.post.mean <- sum.matrices(Z.post) / length(Z.post)
#Z.post.mean <- ifelse(Z.post.mean>runif(length(Z.post.mean)),1,0)
Z.post.mean <- ifelse(Z.post.mean>.5,1,0)
col0.ind <- which(apply(Z.post.mean,2,function(x) sum(x)==0))
Z.post.mean <- Z.post.mean[,-col0.ind]
a.image(Z.post.mean)


one.A <- EAXZ(Y,Z.post.mean,siga=1,sigx=.5)
d2 <- 2
d1 <- ceiling(nrow(one.A)/d2)


pdf("draw.post.out/postA.pdf")
  a.image(one.A,main="Posterior Mean for A")
dev.off()

plot.post.As <- function(one.A) {
  par(mfrow=c(d2,d1),mar=c(.5,.5,1,.5))
  for (i in 1:nrow(one.A)) {
    one.Ai <- matrix(one.A[i,],6,6) # matrix(Y[n,],6,6) = X[[n]]
    a.image(one.Ai,main=paste0("Posterior Mean A",i))
    #a.image(one.Ai,main=paste0("Posterior Mean A",i),col=BLUE)
  }
  par(mfrow=c(1,1))
}

pdf("draw.post.out/postA66.pdf")
  plot.post.As(one.A)
dev.off()

pdf("draw.post.out/postA.pdf")
  par(mfrow=c(1,2))
    a.image(one.A,main="Posterior Mean for A")
    plot.post.As(one.A)
  par(mfrow=c(1,1))
dev.off()

pdf("draw.post.out/Y.pdf")
  a.image(Y,main="X")
dev.off()

pdf("draw.post.out/postZ.pdf")
  a.image(Z.post.mean,main="Posterior Estimate for Z")
dev.off()

pdf("draw.post.out/postAlpha.pdf")
  plot.post(alpha,"Alpha (after 1000 Burn)")
  #plot(density(alpha[-(1:burn)]),main="Posterior for Alpha",col="cornflowerblue",
  #     lwd=3)
dev.off()


post.ZA <- Z.post.mean %*% one.A
plot.post.ZA <- function(n) {
  par(mfrow=c(1,2))
    a.image(matrix(Y[n,],6,6),main=paste0("n=",n,": ",label[n]))
    a.image(matrix(post.ZA[n,],6,6),
            main=paste0("n=",n,":  ",toString(Z.post.mean[n,])))
  par(mfrow=c(1,1))
}

plot.post.each <- function() {
  opts <- par(no.readonly=T)
  par(mfrow=c(5,4),mar=c(.1,.1,1,.1))
    for (n in 10*(1:(nrow(Y)/10))) {
      a.image(matrix(Y[n,],6,6),main=paste0("n=",n,": ",label[n]),cex.main=.8)
      a.image(matrix(post.ZA[n,],6,6),
              main=paste0("n=",n,":  ",toString(Z.post.mean[n,])),cex.main=.8)
     }       
  par(opts)
}

pdf("draw.post.out/postFriends.pdf")
  plot.post.each()
dev.off()

#a.image(matrix(apply(y2,2,mean),6,6))
a.image(Z.post.mean)
plot.post.As(one.A)

A <- matrix(c(0,1,0,0,1,0,
              0,1,0,0,1,0,
              0,0,0,0,0,0,
              1,0,0,0,0,1,
              1,0,0,0,0,1,
              1,1,1,1,1,1),6,6,byrow=T)

pdf("draw.post.out/oneObs.pdf")
  a.image(A)
dev.off()

pdf("draw.post.out/oneObsWithNoise.pdf")
  set.seed(18)
  a.image(A+rnorm(36,0,.5))
dev.off()

pdf("draw.post.out/Z11postpred.pdf")
  X.post.pred <- lapply(as.list((burn+1):length(Z.post)),
                        function(i) Z.post[[i]]%*%EAXZ(Y,Z.post[[i]])+
                                    rnorm(prod(dim(Y)),0,.5))

  X.post.pred.11 <- unlist(lapply(X.post.pred,function(z) z[1,1]))
  plot.post(X.post.pred.11,main="Posterior Predictive for X[1,1]")
dev.off()

pdf("draw.post.out/qq.pdf")
  qq <- mean(X.post.pred.11<Y[1,1])
  #qq <- apply(as.matrix(X.post.pred.11),1,function(x) mean(x<X.post.pred.11))
  plot(density(qq[-1],from=0,to=1),col=col.mult("cornflowerblue","grey80"),lwd=1,
        main="Density of Posterior Predictive Quantiles for X")
  for (i in 1:10) {
    old.time <- Sys.time()
    for (j in 1:ncol(Y)) {
      if (i>1 & j>1) {
        X.pp <- unlist(lapply(X.post.pred,function(z) z[i,j]))
        qq <- apply(as.matrix(X.pp),1,function(x) mean(x<X.pp))
        lines(density(qq[-i],from=0,to=1),col=col.mult("cornflowerblue","grey80"))
        ks.test(qq,"punif")
      }
    }
    count.down(old.time,i,10)
  }
dev.off()

resid <- c(post.ZA-Y)
plot(c(Y),type="l")
lines(c(post.ZA),col="dodgerblue")
lines(resid,col="pink")

pdf("draw.post.out/resid.pdf")
  plot(resid,col="pink",type="l",main="Residuals = PostMean(ZA) - X")
dev.off()
