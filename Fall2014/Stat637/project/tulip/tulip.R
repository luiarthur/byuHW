library(truncnorm) # for rtruncnorm
library(xtable)
library(splines)
source("plotpost.R")

# 12 populations of tulips
#   - 30 measurements made at each of 7 chilling times
#   - Total of 210 = 30*7 measurements per population
#   - none of population 12's germinated (so we exclude)

dat <- read.csv("data/tulip.csv")
dat <- dat[-which(dat$Pop==12),]

y <- ifelse(dat$Germinated=="N",0,1)
N <- length(y)
chill <- dat$ChillingTime
uchill <- unique(chill)
pop <- dat$Population
upop <- unique(pop)
k <- length(unique(pop))
X <- cbind(1,bs(chill)) # bs(chill) => Cubic Splines
colnames(X) <- c("x0","x1","x2","x3")
p <- ncol(X)

create.Z <- function(X,k) {
  p <- ncol(X)
  N <- nrow(X)
  Z <- matrix(0,N,k*p)
  for (i in 1:N) {
    for (j in 1:p) {
      Z[i,pop[i]+k*(j-1)] <- X[i,j]
    }
  }
  Z
}

Z <- create.Z(X,k)

colnames(Z) <- c(paste0("z0.",1:k),paste0("z1.",1:k),paste0("z2.",1:k),paste0("z3.",1:k)) # p=4

system("mkdir -p images")
pdf("images/allPopsData.pdf",width=19,height=13)
    yy <- apply(matrix(uchill),1,function(ct) mean(y[which(chill==ct)]))
    arfplot(uchill,yy,main="All Populations",xlab="Chill Time (Weeks)",
            ylab="Germination Rate",col="grey30",pch=20,cex=4,ylim=c(0,1),
            vgrid=7)
dev.off()


pdf("images/rawData.pdf",width=19,height=13)
  par(mfrow=c(6,2))
    yy <- apply(matrix(uchill),1,function(ct) mean(y[which(chill==ct)]))
    arfplot(uchill,yy,main="All Populations",xlab="Chill Time (Weeks)",
            ylab="Germination Rate",col="grey30",pch=20,cex=2,ylim=c(0,1))
    for (pn in 1:k) {
        yy <- apply(matrix(uchill),1,function(ct) {
          mean(y[which(pop==pn & chill==ct)])
        })
        #plot(uchill,yy,main=paste("Population",pn),ylab="Germination Rate",xlab="Chill Time (Weeks)")
        arfplot(uchill,yy,main=paste("Population",pn),ylab="",xlab="",pch=20,
                col="grey30",cex=2,ylim=c(0,1))
    }
  par(mfrow=c(1,1))
dev.off()


mvrnorm <- function(M,S,n=nrow(S)) M + t(chol(S)) %*% rnorm(n)

updateW <- function(b,g){
  w <- numeric(length(y))
  #w[y==1] <- rtruncnorm(sum(y==1),a=0,b=Inf, mean=X[y==1,]%*%b,sd=1)
  #w[y==0] <- rtruncnorm(sum(y==0),a=-Inf,b=0,mean=X[y==0,]%*%b,sd=1)
  w[y==1] <- rtruncnorm(sum(y==1),a=0,b=Inf, mean=X[y==1,]%*%b + Z[y==1,]%*%g,sd=1)
  w[y==0] <- rtruncnorm(sum(y==0),a=-Inf,b=0,mean=X[y==0,]%*%b + Z[y==0,]%*%g,sd=1)
  w
}

gibbs <- function(y,X,Z,n=nrow(X),p=ncol(X),k=ncol(Z),B=1e4,burn=round(B*.1),
                 trim.burn=F, U=10*solve(t(X)%*%X), V=10*solve(t(Z)%*%Z)) {

  Xt <- t(X)
  Zt <- t(Z)
  XX <- Xt%*%X
  ZZ <- Zt%*%Z
  T1 <- solve(XX + solve(U))
  T2 <- solve(ZZ + solve(V))

  # Initialize Parameters
  b <- matrix(0,B,p)
  g <- matrix(0,B,k)
  colnames(b) <- paste0("b",0:(p-1))
  colnames(g) <- c(paste0("g0.",1:(k/p)),paste0("g1.",1:(k/p)),
                   paste0("g2.",1:(k/p)),paste0("g3.",1:(k/p))) # p=4
  #######################

  for (i in 2:B){
    #Updates:
    old.time <- Sys.time()
    w <- updateW(b[i-1,],g[i-1,])
    #b[i,] <- mvrnorm(T1 %*% Xt %*% w, T1) 
    b[i,] <- mvrnorm(T1 %*% Xt %*% (w-Z%*%g[i-1,]), T1) 

    w <- updateW(b[i,],g[i-1,])
    g[i,] <- mvrnorm(T2 %*% Zt %*% (w-X%*%b[i,]),   T2) 
    count.down(old.time,i,B)
  }

  list("b"=b,"g"=g)
}

out <- gibbs(y,X,Z,B=1e4)
# Plot Posteriors?
#plot.posts(out$b,names=c("b0","b1","b2","b3"))
#plot.posts(out$g[,c(1,11,12)])
#a.image(cov(out$b),axis.num=F)
#a.image(cov(out$g),axis.num=F)
b <- apply(out$b,2,mean)
g <- apply(out$g,2,mean)
B <- out$b
G <- out$g

compare.chill.effect <- function(G) {
  compare.one.pair <- function(i,j) {
    #diff <- G[,i]-G[,j]
    #print(paste(i,j))
    #diff <- apply(G[,c(i+11,i+22,i+33)],1,function(x) sqrt(sum(x^2)))- 
    #        apply(G[,c(j+11,j+22,j+33)],1,function(x) sqrt(sum(x^2)))
    diff <- G[,c(i+11,i+22,i+33)]-G[,c(j+11,j+22,j+33)]
    ci <- t(apply(diff,2,get.hpd))
    ci[1,1] <= 0 && 0 <= ci[1,2] && 
    ci[2,1] <= 0 && 0 <= ci[2,2] &&
    ci[3,1] <= 0 && 0 <= ci[3,2]
  }

  #k <- ncol(G)
  out <- matrix(0,k,k)
  J <- matrix(1:k)
  pairs <- c()
  for (j in 1:k) { # For each population
    ot <- Sys.time()
    for (i in 1:j) {
      out[i,j] <- compare.one.pair(i,j)
      if (out[i,j]==1 && i!=j) pairs <- c(pairs, paste0("(",i,",",j,")"))
    }
    count.down(ot,j,k)
  }

  list("matrix"=out,"pairs"=pairs)
}

comp <- compare.chill.effect(cbind(out$g,out$b[,2]))
# (Last Year) The 9 following pairs behave the most similiarly under different chill times: # (3,2), (4,2), (10,4), (7,6), (10,6), (11,6), (10,7), (11,7), (11,10)
comp$m
comp$p
paste(comp$p,collapse=", ")
# Populations that have essentially the same chilling effects
#(1,2), (2,3), (2,4), (3,4), (4,6), (4,7), (6,7), (1,8), (2,8), (4,8), (5,8), (6,8), (7,8), (6,9), (4,10), (6,10), (7,10), (6,11), (7,11), (10,11)

x0 <- seq(min(uchill),max(uchill),len=100)
X0 <- cbind(1,bs(x0))

post.pred <- function(x,B,G,k=0,f=mean) {
  M <- NULL
  if (k>0) {
    M <- pnorm(x%*%(t(B)+t(G[,c(k,k+11,k+22,k+33)])))
  } else {
    M <- pnorm(x%*%t(B))
  }
  apply(M,1,f)
}

post.pred.all <- function(x,B,G,k=0) {
  M <- NULL
  if (k>0) {
    M <- pnorm(x%*%(t(B)+t(G[,c(k,k+11,k+22,k+33)])))
  } else {
    M <- pnorm(x%*%t(B))
  }
  M
}

best.chill.times <- apply(matrix(1:k),1,function(pn) {
  y0 <- post.pred(X0,B,G,pn)
  pp <- post.pred.all(X0,B,G,pn)
  best.chill.xs <- apply(pp,2,which.max)
  best.chill.x <- x0[mean(best.chill.xs)]
  best.chill.x
})

pop.order <- (1:k)[order(best.chill.times)]


details <- function(pn=0,lwd=1) {
  y0 <- post.pred(X0,B,G,pn)
  pp <- post.pred.all(X0,B,G,pn)
  best.chill.xs <- apply(pp,2,which.max)
  best.chill.x <- x0[mean(best.chill.xs)]
  hpd <- get.hpd(best.chill.xs,len=1e4)
  best.chill.y <- y0[which.min(abs(x0-best.chill.x))]
  points(best.chill.x,best.chill.y,cex=3,col="blue",pch=20)
  #lines(x0,y0,lwd=2,col="grey30")
  lines(x0,y0,lwd=lwd,col="red")
  #lines(x0[which(x0<x0[hpd[1]] | x0>x0[hpd[2]])],
  #      y0[which(x0<x0[hpd[1]] | x0>x0[hpd[2]])],lwd=2,col="grey30")
  lines(x0[which(x0>x0[hpd[1]] & x0<x0[hpd[2]])],
        y0[which(x0>x0[hpd[1]] & x0<x0[hpd[2]])],lwd=lwd,col="blue")
  #lines(c(x0[hpd[1]],x0[hpd[2]]),c(best.chill.y,best.chill.y),col="orange",lwd=3)
  #lines(c(x0[hpd[1]],x0[hpd[1]]),best.chill.y+c(-.1,.1),col="orange",lwd=2)
  #lines(c(x0[hpd[2]],x0[hpd[2]]),best.chill.y+c(-.1,.1),col="orange",lwd=2)
}

pdf("images/chilleffect.pdf",width=19,height=14) # Posterior Predictive Means
  par(mfcol=c(6,2))
    yy <- apply(matrix(uchill),1,function(ct) mean(y[which(chill==ct)])) 
    arfplot(uchill,yy,main="All Populations",xlab="Chill Time (Weeks)",
            ylab="Germination Rate",col="grey30",pch=20,cex=2,ylim=c(0,1),
            vgridlines=7,hgrid=5)
    details(0,lwd=3)
         
    for (pn in pop.order) {
      ot <- Sys.time()
      yy <- apply(matrix(uchill),1,function(ct) {
        mean(y[which(pop==pn & chill==ct)])
      })
      #plot(uchill,yy,main=paste("Population",pn),ylab="Germination Rate",xlab="Chill Time (Weeks)")
      arfplot(uchill,yy,main=paste("Population",pn),ylab="",xlab="",pch=20,
              col="grey30",cex=2,ylim=c(0,1),vgrid=7,hgrid=5)
      details(pn,lwd=3)
      count.down(ot,pn,k)
    }
  par(mfrow=c(1,1))
dev.off()

M <- apply(cbind(out$b,out$g),2,mean)
hpd.param <- t(apply(cbind(out$b,out$g),2,get.hpd))
hpd.param <- hpd.param[order(M),]
M <- M[order(M)]
pdf("images/hpd.pdf",width=19,height=13)
  arfplot(M,ylab="Posterior Estimates",xlab="",mgp=c(1.5,0,0),xaxis=F,vg=48,xlim=range(hpd.param),
          cex.lab=2,cex.yaxis=1.5)
  axis(1,at=1:48,label=names(M),las=2,cex.axis=1.5,tick=F)
  abline(h=0,col="maroon",lwd=2); axis(4,at=0,label=0,las=1,tick=F,cex.axis=1.5)
  apply(matrix(1:48),1,function(x) {
    lines(rep(x,2),hpd.param[x,])
  })
dev.off()

#system("cd latex/report; ./compile")

