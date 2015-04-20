source("rfunctions.R")
source("genData.R")
library(xtable)

system("mkdir -p out")

gibbs.post <- function(y,X,sigg2=100,sigb2=100,a=1,B=1000,burn=B*.1,showProgress=T,a.a=1,a.b=1,a.r=1,b.r=5,cs.r2=2,plotProgress=F) {
  B <- ceiling(B/50)*50 

  D <- ncol(X)
  N <- nrow(X)

  Hn <- sum(1/(1:N))

  p.x.z <- function(Z,sig.b2,sig.g2,sig.r2,log=T) { # p.x.z = Likelihood
    # [y|X,Z,sigb2,sigg2,sigr2]
    #if (ncol(Z) > 10) Z <- Z[,1:10]

    K <- ncol(Z)
    Id <- diag(1,D)
    Ik <- diag(1,K)
    In <- diag(1,N)

    M1 <- solve(t(Z)%*%Z + sig.r2/sig.g2 * Ik)
    U <- In- Z%*%M1%*%t(Z)
    u <- chol(U)
    ys <- u%*%y
    xs <- u%*%X
    M2 <- solve(t(xs)%*%xs + sig.r2/sig.b2 * Id)

    out <- 0
    if (!log) {
      out <- (2*pi*sig.r2)^(-N/2) * (sig.r2/sig.g2)^(K/2) * (sig.r2/sig.b2)^(D/2) *
             det(M1) * det(M2) * 
             exp(-.5/sig.r2 * t(ys) %*% (In-xs%*%M2%*%t(xs)) %*% ys)
    } else {
      out <- -N/2*log(2*pi*sig.r2)+K/2*log(sig.r2/sig.g2)+D/2*log(sig.r2/sig.b2)+
             det(M1,log=T)+det(M2,log=T) +
             -.5/sig.r2 * t(ys) %*% (In-xs%*%M2%*%t(xs)) %*% ys
    }

    out
  } # p.x.z
  
  p.zik.x <- function(z,i,k,sig.b2,sig.g2,sig.r2,log=T) { # P[z_{ik}=1|z_{-i,k}].  
                                                          # exact

    zz <- z
    zz[i,k] <- 1
    a <- p.x.z(zz,sig.b2,sig.g2,sig.r2,log=T) 
    zz[i,k] <- 0
    b <- p.x.z(zz,sig.b2,sig.g2,sig.r2,log=T) 
    mk <- sum(z[-i,k])
    p1 <- mk/N # p = Prior
    p0 <- 1 - p1

    p0 <- b+log(p0) 
    p1 <- a+log(p1)
    
    if (mk==0){
      out <- 0 # essentially doing nothing. because mk=0 => z[i,k] = 0.
    } else if (!log) { # Not Logged
      out <- 1 / (1+exp(p0-p1))
    } else { # Logged
      out <- -log(1+exp(p0-p1))
    }
    
    out
  }
  

  sampNewCols <- function(z,i,s=9,a=alpha[1],sig.b2,sig.g2,sig.r2) { #s is the max num of columns. Avoid mh.
    prior <- function(l) dpois(l,a/N,log=T)  

    lp <- prior(0:s) # log prior
    ll <- apply(matrix(0:s),1,function(x) { 
                                col1 <- matrix(0,N,x)
                                col1[i,] <- 1
                                p.x.z(cbind(z,col1),sig.b2,sig.g2,sig.r2,log=TRUE)
                                }) # log like
    lg <- lp+ll

    lpi <- apply(matrix(0:s),1,function(i) -log(sum(exp(lg-lg[i+1]))))
    Pi <- exp(lpi)

    out <- sample(0:s,1,prob=Pi)
    out
  }

  Zs <- as.list(1:B)
  Zs[[1]] <- matrix(1,N,1)
  alpha <- rep(a,B)
  sig.r2 <- rep(1,B)
  count.acc.r2 <- 0

  for (b in 2:B) { # B = num of iterations in Gibbs
    old.time <- Sys.time()
    z <- Zs[[b-1]]
    alpha[b] <- alpha[b-1]
    sig.r2[b] <- sigr2 <- sig.r2[b-1]

    for (i in 1:N) { # iterate through all rows of Z
      cat("\ri=",i)
      K <- ncol(z)

      k <- 1
      while(k<=K) { # iterate through all columns of Z
        if (K>0) {
          p <- p.zik.x(z,i,k,sigb2,sigg2,sigr2,log=F) # pzik=1|x.
          u <- runif(1)
          if (p<=0) {
            z <- as.matrix(z[,-k])
            k <- k-1
          } else if (p>u){
            z[i,k] <- 1
          } else {
            z[i,k] <- 0
          }
        }

        k <- k + 1
        K <- ncol(z)
      } # end of its through all columns of Z

      # Drawing new dishes should be prior times likelihood
      new <- sampNewCols(z,i,s=15,a=alpha[b],sigb2,sigg2,sigr2) # Original: March 14, 2015
      #new <- rpois(1,alpha[b]/N) # New: March 14, 2015

      if (new>0) {
        col1 <- matrix(0,N,new)
        col1[i,] <- 1
        z <- as.matrix(cbind(z,col1))
      }
    } # end of its through all rows of Z
# Z|alpha propto alpha^{a-1} exp{-alpha Hn} #   alpha propto alpha^{a-1} exp{-alpha / b}
    # if a=1, b=1
    # a|Z ~ Gamma(a+K,(1/b+Hn)^(-1)) 
    #a.a <- a.a+ncol(z)
    #a.b <- (1/a.b+Hn)^(-1)
    alpha[b] <- rgamma(1,a.a+ncol(z),scale=1/(1/a.b+Hn))
    Zs[[b]] <- z

    # MH: Sample sig.r2. VERSION 3:
    cand <- rnorm(1,sigr2,cs.r2)
    if (cand>0) {
      lr2.cand <- -(a.r-1)*log(cand)-b.r/cand+p.x.z(z,sigb2,sigg2,cand,log=T)
      lr2.old <- -(a.r-1)*log(sigr2)-b.r/sigr2+p.x.z(z,sigb2,sigg2,sigr2,log=T)
      lr2 <- lr2.cand-lr2.old
      if (lr2 > log(runif(1))) {
        sig.r2[b] <- cand
        count.acc.r2 <- count.acc.r2 + 1
      }
    }

    # Print Results as I go
    sink("out/Z.post.results",append=b>2)
      cat("ITERATION:",b,"\n")
      print(z)
    sink()

    if (plotProgress && b%%10==0) {
      # Plot Trace Plots:
      n.col <- unlist(lapply(Zs[1:b],ncol))
      plot(n.col,xlab="Iteration",ylab="K+",
           main=paste0("Columns of Z ","(",b,")"),col="pink",lwd=3,type="b",pch=20)
      abline(h=mean(n.col),col="blue",lwd=3)     
      minor <- function() {plot(alpha[1:b],type="l",col="gray30",cex.main=.6,
                                main=paste("alpha:",round(mean(alpha[1:b]),4))) } 
      minor2 <- function() {plot(sig.r2[1:b],type="l",col="gray30",cex.main=.6,
                                 main=paste("sig.r2:",round(mean(sig.r2[1:b]),4)))}
      plot.in.plot(minor,"topright")
      plot.in.plot(minor2,"bottomright")

    }

    if (showProgress) count.down(old.time,b,B)
  } # end of gibbs

  cat("\n Acceptance for sig.r2: ",count.acc.r2/B,"\n")
  #out <- Zs[(burn+1):B]
  out <- list("Zs"=Zs,"alpha"=alpha,"sig.r2"=sig.r2)
  out
}

# SIMULATIONS STUDY: UNCOMMENT TO SIMULATE!!!
# Works when sig.r2=1
  dat <- gendata()
  y <- dat$y
  X <- dat$X
  b <- dat$b
  Z <- dat$Z
  g <- dat$g
  sigr2 <- dat$sigr2
# End of genData.R Plots

# Simulation:
B <- 1e3
elapsed.time <- system.time(out <- gibbs.post(y,X,cs.r2=.2,B=B,showProgress=T,plotProgress=F))

# Analysis:
system("mkdir -p latex/images")

EZ.pre <- est.Z(out$Zs)
EZ <- clust.Z(EZ.pre)

pdf("latex/images/EZpre.pdf")
  a.image(EZ.pre,axis.num=F,col=paste0("grey",100:30))
dev.off()

 
pdf("latex/images/EZ.pdf")
  a.image(EZ,axis.num=F,color=paste0("grey",100:30))
dev.off()

pdf("latex/images/posta.pdf")
  plot.post(out$a)
dev.off()

pdf("latex/images/postsigr2.pdf")
  plot.post(out$sig.r2)
dev.off()

pdf("latex/images/postMeanZ.pdf")
  a.image(sum.matrices(out$Zs)/B,axis.num=F,color=paste0("grey",100:30))
dev.off()

pdf("latex/images/scatter.pdf")
  plot(X[,2],y,pch=20,xlab="x")
dev.off()  

pdf("latex/images/clus.pdf")
  plot.mm(y,X[,2],b,Z,g,pch=20,line=F)
dev.off()  

pdf("latex/images/agpost.pdf")
  plot.posts(cbind(out$a,out$sig.r2),names=c("a","sig.r2"))
dev.off()

G <- diag(100,ncol(EZ))
R <- diag(mean(out$sig.r2),length(y))
V <- EZ %*% G %*% t(EZ) + R
beta.hat <- solve(t(X) %*%solve(V) %*%X)%*%t(X) %*%solve(V) %*%y
gam.hat <- G%*%t(EZ)%*%solve(V)%*%(y-X%*%beta.hat)

pdf("latex/images/resultmm.pdf")
  plot.mm(y,X[,2],beta.hat,EZ,gam.hat,pch=20)
dev.off()

pdf("latex/images/truemm.pdf")
  plot.mm(y,X[,2],b,EZ,g,pch=20)
dev.off()

kmean.clus <- kmeans(cbind(y,X[,2]),3)$cluster+4
kmean.clus <- ifelse(kmean.clus==6,4,ifelse(kmean.clus==5,2,3))
pdf("latex/images/kmean.pdf")
  plot(X[,2],y,xlab="x",col=kmean.clus,pch=20)
dev.off()

KZ <- matrix(0,n,3)
for (i in 1:n) {
  KZ[i,kmean.clus[i]-1] <- 1
}
KZ


Gk <- diag(100,ncol(KZ))
Rk <- diag(mean(out$sig.r2),length(y))
Vk <- KZ %*% Gk %*% t(KZ) + R
bk.hat <- solve(t(X) %*%solve(Vk) %*%X)%*%t(X) %*%solve(Vk) %*%y
gk.hat <- Gk%*%t(KZ)%*%solve(Vk)%*%(y-X%*%bk.hat)

pdf("latex/images/ibpmm.pdf")
  plot.mm(y,X[,2],b=beta.hat,z=EZ,gam.hat,pch=20)
dev.off()

pdf("latex/images/KMmm.pdf")
  plot.mm(y,X[,2],b=bk.hat,z=KZ,gk.hat,pch=20)
dev.off()

sink("latex/images/mb.tex")
  Mb <- cbind(b,beta.hat,bk.hat)
  rownames(Mb) <- c("$\\beta_0$","$\\beta_1$")
  colnames(Mb) <- c("True","IBP","KM")
  Mb <- xtable(Mb,digits=c(0,0,3,3))
  print(Mb,sanitize.text.function=function(x) x)
sink()

sink("latex/images/mg.tex")
  Mg <- cbind(g,gam.hat,gk.hat)
  rownames(Mg) <- c("$\\gamma_0$","$\\gamma_1$","$\\gamma_2$")
  colnames(Mg) <- c("True","IBP","KM")
  Mg <- xtable(Mg,digits=c(0,0,3,3))
  print(Mg,sanitize.text.function=function(x) x)
sink()

pdf("combineMM.pdf",height=11)
  par(mfrow=c(3,1))
    plot.mm(y,X[,2],b=b,z=Z,g,pch=20,main="TRUE")
    plot.mm(y,X[,2],b=beta.hat,z=EZ,gam.hat,pch=20,main="IBP")
    plot.mm(y,X[,2],b=bk.hat,z=KZ,gk.hat,pch=20,main="K-Means")
  par(mfrow=c(1,1))
dev.off()

pdf("latex/images/KZ.pdf")
  a.image(KZ,axis.num=F,col=paste0("grey",100:30))
dev.off()
 
