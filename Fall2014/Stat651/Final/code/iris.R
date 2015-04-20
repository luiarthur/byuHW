source("rfunctions.R")
source("ibp.R")
source("gibbs.R")

Y <- as.matrix(iris[,1:4])
label <- as.vector(iris[,5])
elapsed.time <- system.time(out <- gibbs.post(Y,a=1,B=1000,burn=0,showProgress=T,
                                              plotProgress=T,a.a=1,a.b=2,
                                              siga=1,sigx=1))

M <- out$Zs
alpha <- out$alpha
#burn <- 1000#round(length(M) * .1)
burn <- round(length(M) * .1)


n.col <- unlist(lapply(M,ncol))

plot.trace.ncol <- function() {  
  plot(n.col,type="l",main="Trace Plot: Number of Columns in Z",lwd=1,cex=.1,
       col="gray30",pch=20)
  mean.col <- round(mean(n.col[-(1:burn)]),4)
  var.col <-  round(var(n.col[-(1:burn)]),5)
  abline(h=mean.col,lwd=2,col="red")
  legend("bottomright",legend=c(paste("Mean=",mean.col), 
                                paste("Variance =",var.col)),title.col="gray30",
                                title=paste("After Burn-in of",burn,":"),bty="n")
}                                

pdf("out/traceplot.pdf")
  plot.trace.ncol()
dev.off()


plot.post.alpha <- function() {  
  plot(alpha,type="l",main="Trace Plot: Alpha",lwd=1,cex=.1,
       col="gray30",pch=20)
  mean.a <- round(mean(alpha[-(1:burn)]),4)
  var.a <-  round(var(alpha[-(1:burn)]),5)
  abline(h=mean.a,lwd=2,col="red")
  legend("topleft",legend=c(paste("Mean=",mean.a),
                                paste("Variance =",var.a)),title.col="gray30",
                                title=paste("After Burn-in of",burn,":"),bty="n")
}                                

pdf("out/tracealpha.pdf")
  plot.post.alpha()
dev.off()

EAXZ <- function(X,Z,siga=1,sigx=1) {
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


