source("countdown.R")
source("../3/color.R")
library(xtable)
options("width"=100)

system("mkdir -p out")
y <- as.vector(as.matrix(read.table("faculty.dat")))
k <- length(y)

dig <- function(x,a,b,log=F) { # b is rate. 1/b is scale.
  out <- NULL
  if (!log) {
    out <- 1/(gamma(a)*b^a) * x^(-a-1) * exp(-1/(b*x))
  } else {
    out <- -lgamma(a)-a*log(b) -(a+1)*log(x) -1/(b*x)
  }

  out
}

rig <- function(n,a,b) 1/rgamma(n,a,scale=b)

# Prior Speicifications
  n <- 1e5
  #Priors:
  m <- 4.5   # => E[mu] = 4.5
  s2 <- 2    # => V[mu] = 2
  as <- 5/2  # => E[sig2] = 1
  bs <- 2/3  # => V[sig2] = 2
  at <- 13/2 # => E[tau2] = 3
  bt <- 2/33 # => V[tau2] = 2

  prior.tau2 <- rig(n,at,sqrt(bt))
  prior.sig2 <- rig(n,as,sqrt(bs))
  prior.mu <- rnorm(n,m,sqrt(s2))
  prior.theta <- rnorm(n,prior.mu,sqrt(prior.tau2))
  prior.y <- rnorm(n,prior.theta,sqrt(prior.sig2))

  pdf("out/priorPred.pdf")
    plot(density(prior.y),lwd=3,col="red",main="Prior Predictive")
    color.den(density(prior.y),1,7,"red")
    #text(6,.1,round(mean(prior.y>5),4),cex=2)
    text(4.5,.1,round(mean(prior.y>1 & prior.y<7),4),cex=2)
    #abline(v=mean(prior.y)) # mean about 4.5, var about 3.56
  dev.off()

gibbs <- function(B=1e5) {
  M <- matrix(0,nrow=B,ncol=k+3) # theta[1:k],mu,sig2,tau2
  M[1,] <- 2
  
  update.theta <- function(i,mu,sig2,tau2) { # theta ~ N(mu,sig2)
    mu.new <- (y[i]*tau2+mu*sig2) / (tau2+sig2)
    sig2.new <- (tau2*sig2) / (tau2+sig2)
    out <- c(mu.new,sig2.new)
    names(out) <- c("mu","sig2")
    out
  }

  update.mu <- function(theta,tau2) { # mu ~ N(m,s2)
    theta.new <- (mean(theta)*k*s2+m*tau2)/(k*s2+tau2)
    s2.new <- s2*tau2/(k*s2+tau2)
    out <- c(theta.new,s2.new)
    names(out) <- c("m","s2")
    out
  }

  update.sig2 <- function(theta) { # sig2 ~ IG(as,bs)
    as.new <- as+k/2 # Does not change
    bs.new <- 1/(1/bs+sum((y-theta)^2)/2)
    out <- c(as.new,bs.new)
    names(out) <- c("as","bs")
    out
  }
  
  update.tau2 <- function(theta,mu) { # tau2 ~ IG(at,bt)
    at.new <- at+k/2 # Does not change
    bt.new <- 1/(1/bt+sum((theta-mu)^2)/2)
    out <- c(at.new,bt.new)
    names(out) <- c("at","bt")
    out
  }

  for (i in 2:B) {
    old.time <- Sys.time()
    M[i,] <- M[i-1,]
    
    # theta[1:k],mu,sig2,tau2
    for (j in 1:k) {
      new.theta <- update.theta(j,M[i,k+1],M[i,k+2],M[i,k+3])
      M[i,j] <- rnorm(1,new.theta[1],new.theta[2])
    }

    new.mu <- update.mu(M[i,1:k],M[i,k+3])
    M[i,k+1] <- rnorm(1,new.mu[1],new.mu[2])

    new.sig2 <- update.sig2(M[i,1:k])
    M[i,k+2] <- rig(1,new.sig2[1],new.sig2[2])

    new.tau2 <- update.tau2(M[i,1:k],M[i,k+1])
    M[i,k+3] <- rig(1,new.tau2[1],new.tau2[2])
    
    #print(M[i,(k+1):(k+3)])
    #print(new.tau2)
    #Sys.sleep(1)
    #cat(paste0("\r",round(i/B*100),"%"))
    count.down(old.time,i,B)
  }
  
  M
}

#1: Posterior
M <- gibbs(B=1e5)
M <- M[-(1:200),]
N <- nrow(M)

#for (i in 1:ncol(M)){
#  plot(M[,i],main=i,type="l")
#  scan()
#}


#2: E[theta|Y]
M.mean <- apply(M,2,mean)
M.upper <- apply(M,2,function(x) quantile(x,.975))
M.lower <- apply(M,2,function(x) quantile(x,.025))

#3: V[theta|Y]
M.cov <- var(M)
#sum(M.cov > .05) / prod(dim(M.cov))

# Plots for mean and variance
pdf("out/postMean.pdf")  
  plot(M.mean[1:23],col="blue",pch=20,
       ylim=c(min(M.lower[1:23]),max(M.upper[1:23])),
       type="o",main=expression(paste("Posterior Mean ",theta,"'s")),xlab="i",
       ylab=expression(theta[i]))
  lines(M.upper[1:23],col="red",type="b",pch=20)
  lines(M.lower[1:23],col="red",type="b",pch=20)
dev.off()

sink("out/postMeanMu.tex")
  P <- cbind(M.mean[24:26],M.lower[24:26],M.upper[24:26],diag(M.cov)[24:26])
  colnames(P) <- c("Mean","CI.lower","CI.upper","Variance")
  rownames(P) <- c("$\\mu$","$\\sigma^2$","$\\tau^2$")
  print(xtable(P,digits=6,caption="Posterior Mean \\& Variance"),sanitize.t=function(x)x,caption.p="top")
sink()

library(MASS)
plot.hyper <- function(m,names=NULL) {
  k <- ncol(m)
  par(mfrow=c(k,k))
  for (i in 1:k) {
    for (j in 1:k) {
      if (i<j) {
        name <- ifelse(is.null(names),"",paste("Bivariate Contour Plot for \n",names[i],"&",names[j]))
        K <- kde2d(m[,i],m[,j])
        contour(K,col="red",main=name,cex.main=.9,xlab=names[i],ylab=names[j])
      } else if (i>j) {
        name <- ifelse(is.null(names),"",
                paste("Trace Plot Between \n",names[i],"&",names[j]))
        plot(m[,i],m[,j],type="l",col="pink",main=name,cex.main=.9,
             xlab=names[i],ylab=names[j])
        points(m[1,i],m[1,j],pch="B")     
        points(m[nrow(m),i],m[nrow(m),j],pch="E")     
        #plot.new()
        #par(usr=c(0,1,0,1))
        #text(.5,.5,round(cov(m[,i],m[,j]),4))
        #box()
        #name <- ifelse(is.null(names),"",
        #        paste("Covariance between \n",names[i],"&",names[j]))
        #title(name,cex.main=.9)

      } else {
        name <- ifelse(is.null(names),"",paste("Posterior Density for",names[i]))
        plot(density(m[,i]),col="blue",lwd=3,main=name,cex.main=.9,xlab=names[i])
      }
    }
  }
  par(mfrow=c(1,1))
}



pdf("out/postDist.pdf")
  plot.hyper(m<-M[(nrow(M)*.8):nrow(M),24:26],c("mu","sigma2","tau2"))
  #plot.hyper(m<-M[(nrow(M)*.8):nrow(M),22:26],c("theta22","theta23","mu","sigma2","tau2"))
  #plot.hyper(m<-M[(nrow(M)*.8):nrow(M),23:26],c("theta23","mu","sigma2","tau2"))
dev.off()

sink("out/meanVar.tex")
  mv <- cbind(M.mean,diag(M.cov))
  rownames(mv) <- c(paste0("$\\theta_{",1:23,"}$"),
                    "$\\mu$","$\\sigma^2$","$\\tau^2$")
  colnames(mv) <- c("Mean","Variance")
  print(xtable(mv,digits=6,caption="Posterior Mean \\& Variance"),
        sanitize.t=function(x)x,caption.p="top")
sink()

pdf("out/postDistUni.pdf")
  plot(density(M[,1]),col="blue",lwd=3,
       expression(paste("Posterior Distribution For ",theta[1])))
dev.off()  
pdf("out/traceUni.pdf")  
  plot(M[,1],type="l",col="purple",
       main=expression(paste("Trace Plot For ",theta[1])))
dev.off()

pdf("out/trace.pdf")
  par(mfrow=c(4,1))
    plot(M[,23],type="l",col="purple",ylab=expression(theta[23]),
         main=expression(paste("Trace Plot For ",theta[23])))
    plot(M[,24],type="l",col="purple",ylab=expression(mu),
         main=expression(paste("Trace Plot For ",mu)))
    plot(M[,25],type="l",col="purple",ylab=expression(sigma^2),
         main=expression(paste("Trace Plot For ",sigma^2)))
    plot(M[,26],type="l",col="purple",ylab=expression(tau^2),
         main=expression(paste("Trace Plot For ",tau^2)))
  par(mfrow=c(1,1))
dev.off()

#4: Posterior Predictive
#                     mu           tau2
theta.pred <- rnorm(N,M[,k+1],sqrt(M[,k+3]))
#                                    sig2
post.pred <- rnorm(N,theta.pred,sqrt(M[,k+2]))
#post.pred <- ifelse(post.pred>7,7,post.pred)
#post.pred <- ifelse(post.pred<0,0,post.pred)
#post.pred.den <- density(post.pred,from=0,to=7)
post.pred.den <- density(post.pred)
p.gt.5 <- mean(post.pred>5)
mx <- max(post.pred.den$x)

pdf("out/postPred.pdf")
  plot(post.pred.den,lwd=3,col="blue",main="Posterior Predictive for Next Average Faculty Evaluation")
  color.den(post.pred.den,5,mx,col="blue")
  text(5.8,.2,round(p.gt.5,4),cex=2)
dev.off()

a.image <- function(Q,color=paste0("gray",0:100),...) {
  library(fields)
  pr <- par("mar")
  par(mar=c(5,4.5,4,7))
  image(t(apply(Q,2,rev)),yaxt="n",xaxt="n",col=color,...)
  axis(2,at=seq(0,1,length.out=ncol(Q)),labels=rev(colnames(Q)),las=2)
  axis(3,at=seq(0,1,length.out=nrow(Q)),labels=rownames(Q),las=2)
  image.plot(Q,legend.only=T,col=color)
  par(mar=pr)
}

cov.M <- cov(M)
rownames(cov.M) <- c(paste0("theta",1:23),"mu","sig2","tau2")
colnames(cov.M) <- rownames(cov.M)

pdf("out/cov.pdf")
  a.image(cov.M,col=tim.colors(12))
dev.off()


