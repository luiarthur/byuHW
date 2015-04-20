source("rfunctions.R")

gibbs.post <- function(y,X,siga=1,sigx=.5,a=1,B=1000,burn=B*.1,showProgress=T,
                       a.a=1,a.b=1,sg.a=1,sg.b=1,sr.a=1,sr.b=1,plotProgress=F) {
  B <- ceiling(B/50)*50 

  D <- ncol(X)
  N <- nrow(X)

  Hn <- sum(1/(1:N))

  p.x.z <- function(Z,bb,sig.g2,sig.r2,log=T) { # p.x.z = Likelihood
    K <- ncol(Z)
    G <- diag(sig.g2,K)
    R <- diag(sig.r2,N)
    V <- Z %*% G %*% t(Z) + R # V = ZGZ'+R

    A <- y-X%*%bb
    if (!log) {
      out <- (2*pi)^(-N/2)*exp(-.5*t(A) %*% solve(V) %*% A) / sqrt(det(V))
    } else {
      out <- -N/2*log(2*pi) -.5 * (t(A)%*% solve(V) %*% A + log(det(V)))
    }
    out
  } # p.x.z
  
  p.zik.x <- function(z,i,k,b,sig.g2,sig.r2,log=T) { # P[z_{ik}=1|z_{-i,k}].  # exact

    zz <- z
    zz[i,k] <- 1
    a <- p.x.z(zz,b,sig.g2,sig.r2,log=T) 
    zz[i,k] <- 0
    b <- p.x.z(zz,b,sig.g2,sig.r2,log=T) 
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
  

  sampNewCols <- function(z,i,s=9,a=alpha[1],b.hat,sig.g2,sig.r2) { #s is the max num of columns. Avoid mh.
    prior <- function(l) dpois(l,a/N,log=T)  

    lp <- prior(0:s) # log prior
    ll <- apply(matrix(0:s),1,function(x) { 
                                col1 <- matrix(0,N,x)
                                col1[i,] <- 1
                                p.x.z(cbind(z,col1),b.hat,sig.g2,sig.r2)
                                }) # log like
    lg <- lp+ll

    lpi <- apply(matrix(0:s),1,function(i) -log(sum(exp(lg-lg[i+1]))))
    Pi <- exp(lpi)

    out <- sample(0:s,1,prob=Pi)
    out
  }

  Zs <- as.list(1:B)
  Zs[[1]] <- matrix(1,N,0)
  alpha <- rep(a,B)
  sig.g2 <- rep(1,B)
  sig.r2 <- rep(1,B)
  b.hat <- matrix(0,B,D)
  gam <- as.list(1:B)

  for (b in 2:B) { # B = num of iterations in Gibbs
    old.time <- Sys.time()
    z <- Zs[[b-1]]
    alpha[b] <- alpha[b-1]
    b.hat[b,] <- bhat <- b.hat[b-1,]
    gam[[b]] <- gam[[b-1]]
    sig.g2[b] <- sigg2 <- sig.g2[b-1]
    sig.r2[b] <- sigr2 <- sig.r2[b-1]

    for (i in 1:N) { # iterate through all rows of Z
      cat("\ri=",i)
      K <- ncol(z)

      k <- 1
      while(k<=K) { # iterate through all columns of Z
        if (K>0) {
          p <- p.zik.x(z,i,k,bhat,sigg2,sigr2,log=F) # pzik=1|x.
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
      new <- sampNewCols(z,i,s=15,a=alpha[b],bhat,sigg2,sigr2) # Original: March 14, 2015
      #new <- rpois(1,alpha[b]/N) # New: March 14, 2015

      if (new>0) {
        col1 <- matrix(0,N,new)
        col1[i,] <- 1
        z <- as.matrix(cbind(z,col1))
      }
    } # end of its through all rows of Z

    # Z|alpha propto alpha^{a-1} exp{-alpha Hn}
    #   alpha propto alpha^{a-1} exp{-alpha / b}
    # if a=1, b=1
    # a|Z ~ Gamma(a+K,(1/b+Hn)^(-1)) 
    #a.a <- a.a+ncol(z)
    #a.b <- (1/a.b+Hn)^(-1)
    alpha[b] <- rgamma(1,a.a+ncol(z),scale=1/(1/a.b+Hn))
    Zs[[b]] <- z
    sig.g2[b] <- 73 # should draw from posterior instead of set to 1 every time!!!
    sig.r2[b] <- 2 # should draw from posterior instead of set to 1 every time!!!
    G <- diag(sig.g2[b],ncol(z))
    R <- diag(sig.r2[b],N)
    V <- z %*% G %*% t(z) + R
    Vinv <- solve(V)
    tXVinv <- t(X) %*% Vinv
    b.hat[b,] <- solve(tXVinv %*% X) %*% tXVinv %*% y
    gam[[b]] <- G %*% t(z) %*% Vinv %*% (y-X%*%b.hat[b,])
    
    sink("out/Z.post.results",append=b>2)
      cat("ITERATION:",b,"\n")
      print(b.hat[b,])
      print(z)
      print(gam[[b]])
    sink()
    #if (b %% 50 == 0) {
    #  gp <- b%/%50
    #  sink("out/Z.post.results",append=T)
    #    Zs[ ((gp-1)*50+1) : (gp*50) ]
    #  sink()  
    #} else if (b==1) {
    #  sink("out/Z.post.results")
    #  sink()  
    #}

    if (plotProgress && b%%10==0) {
      # Plot Trace Plots:
      n.col <- unlist(lapply(Zs[1:b],ncol))
      plot(n.col,xlab="Iteration",ylab="K+",
           main=paste0("Columns of Z ","(",b,")"),col="pink",lwd=3,type="b",pch=20)
      abline(h=mean(n.col),col="blue",lwd=3)     
      minor <- function() { plot(alpha[1:b],type="l",col="gray30",main="Trace for alpha") } 
      plot.in.plot(minor,"topright")

    }

    if (showProgress) count.down(old.time,b,B)
  } # end of gibbs

  #out <- Zs[(burn+1):B]
  out <- list("Zs"=Zs,"alpha"=alpha,"gam"=gam,"b"=b.hat,"sig.g2"=sig.g2,"sig.r2"=sig.r2)
  out
}

# SIMULATION: UNCOMMENT TO SIMULATE!
#source("genData.R")
#
#B <- 10000
#elapsed.time <- system.time(out <- gibbs.post(y,X,B=B,showProgress=T,plotProgress=T))
#
#EZ <- est.Z(out$Zs)
##EZ <- Z
#a.image(EZ,axis.num=F,main="Posterior Estimate of Z")
#
#G <- diag(mean(out$sig.g2),ncol(EZ))
#R <- diag(mean(out$sig.r2),length(y))
#V <- EZ %*% G %*% t(EZ) + R
#beta.hat <- solve(t(X) %*%solve(V) %*%X)%*%t(X) %*%solve(V) %*%y
#gam.hat <- G%*%t(EZ)%*%solve(V)%*%(y-X%*%beta.hat)
#
#plot(X,y)
#points(X,X%*%beta.hat+EZ%*%gam.hat,col="blue",cex=2)
#
#cbind(b,beta.hat)
#gam
#gam.hat
