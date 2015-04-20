source("rfunctions.R")
source("ibp.R")
source("gibbs.R")

#Y <- as.matrix(read.table("phil3.dat"))
source("genDat.R")
Y <- Y[1:100,]
elapsed.time <- system.time(out <- gibbs.post(Y,a=1,B=2000,burn=0,showProgress=T,
                                              plotProgress=T,a.a=3,a.b=2,
                                              siga=1,sigx=.5))


M <- out$Zs
alpha <- out$alpha
burn <- 1000#round(length(M) * .1)

n.col <- unlist(lapply(M,ncol))

pdf("out/traceplot.pdf")
  plot(n.col,type="l",main="Trace Plot: Number of Columns in Z",lwd=1,cex=.1,
       col="gray30",pch=20)
  mean.col <- round(mean(n.col[-(1:burn)]),4)
  var.col <-  round(var(n.col[-(1:burn)]),5)
  abline(h=mean.col,lwd=2,col="red")
  legend("bottomright",legend=c(paste("Mean=",mean.col), 
                                paste("Variance =",var.col)),title.col="gray30",
                                title=paste("After Burn-in of",burn,":"),bty="n")
dev.off()

pdf("out/tracealpha.pdf")
  plot(alpha,type="l",main="Trace Plot: Alpha",lwd=1,cex=.1,
       col="gray30",pch=20)
  mean.a <- round(mean(alpha[-(1:burn)]),4)
  var.a <-  round(var(alpha[-(1:burn)]),5)
  abline(h=mean.a,lwd=2,col="red")
  legend("topleft",legend=c(paste("Mean=",mean.a),
                                paste("Variance =",var.a)),title.col="gray30",
                                title=paste("After Burn-in of",burn,":"),bty="n")
dev.off()

EAXZ <- function(X,Z,siga=1,sigx=sigX) {
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
if (length(col0.ind)>0) Z.post.mean <- Z.post.mean[,-col0.ind]
a.image(Z.post.mean)


one.A <- EAXZ(Y,Z.post.mean,siga=1,sigx=.5)
d2 <- 2
d1 <- ceiling(nrow(one.A)/d2)


pdf("out/postA.pdf")
  a.image(one.A,main="Posterior Mean for A")
dev.off()

plot.post.As <- function(one.A) {
  par(mfrow=c(d2,d1),mar=c(.5,.5,1,.5))
  for (i in 1:nrow(one.A)) {
    one.Ai <- matrix(one.A[i,],6,6) # matrix(Y[n,],6,6) = X[[n]]
    a.image(one.Ai,main=paste0("Latent Feat.",i),cex.main=.8)
    #a.image(one.Ai,main=paste0("Posterior Mean A",i),col=BLUE)
  }
  par(mfrow=c(1,1))
}

pdf("out/postA66.pdf")
  plot.post.As(one.A)
dev.off()

pdf("out/Y.pdf")
  a.image(Y,main="X")
dev.off()

pdf("out/postZ.pdf")
  a.image(Z.post.mean,main="Posterior Estimate for Z")
dev.off()

pdf("out/postAlpha.pdf")
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
      a.image(matrix(Y[n,],6,6),main=paste0("n=",n,": ",toString(Z[n,])),cex.main=.8)
      a.image(matrix(post.ZA[n,],6,6),
              main=paste0("n=",n,":  ",toString(Z.post.mean[n,])),cex.main=.8)
     }       
  par(opts)
}

pdf("out/postFriends.pdf")
  plot.post.each()
dev.off()


#a.image(matrix(apply(y2,2,mean),6,6))
a.image(Z.post.mean)
plot.post.As(one.A)
a.image(matrix(post.ZA[70,],6,6))


pdf("out/Z11postpred.pdf")
  X.post.pred <- lapply(as.list((burn+1):length(Z.post)),
                        function(i) Z.post[[i]]%*%EAXZ(Y,Z.post[[i]])+
                                    rnorm(prod(dim(Y)),0,.5))
  X.post.pred.11 <- unlist(lapply(X.post.pred,function(z) z[1,1]))
  plot.post(X.post.pred.11,main="Posterior Predictive for Z[1,1]")
dev.off()

pdf("out/qq.pdf")
  qq <- apply(as.matrix(X.post.pred.11),1,function(x) mean(x<X.post.pred.11))
  plot(density(qq,from=0,to=1),col=col.mult("cornflowerblue","grey80"),lwd=5,
       main="Density of Posterior Predictive Quantiles for Z")
  ks <- ks.test(qq,"punif")
  legend("bottomleft",legend=paste("P-val for KS-stat =",ks$p.value),bty="n")
dev.off()     

