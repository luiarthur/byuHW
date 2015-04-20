source("rfunctions.R")

gibbs.post <- function(X=Y,siga=1,sigx=.5,a=1,B=1000,burn=B*.1,showProgress=T,
                       a.a=1,a.b=1,plotProgress=F) {
  B <- ceiling(B/50)*50 

  D <- ncol(X)
  N <- nrow(X)

  Z <- rIBP(N,.5)
  alpha <- NULL
  alpha[1] <- a

  Hn <- sum(1/(1:N))

  p.x.z <- function(Z,log=T,with.norm.const=T) { # p.x.z = Likelihood
    K <- ncol(Z)
    H <- Z %*% solve(t(Z)%*%Z+(sigx/siga)^2 * diag(K)) %*% t(Z)
    I <- diag(nrow(H))

    if (!log) {
      norm.const <- ifelse(with.norm.const,
                           (2*pi)^(N*D/2) * sigx^((N-K)*D) * siga^(K*D),1)
      out <- 
      ( norm.const * det(t(Z)%*%Z + (sigx/siga)^2 * diag(K))^(D/2) )^(-1) *
      exp(
        -1/(2*sigx^2) * tr(t(X) %*% (I-H) %*% X)
      )
    } else {
      norm.const <- ifelse(with.norm.const,
                           (N*D/2)*log(2*pi)+((N-K)*D)*log(sigx)+(K*D)*log(siga),0)
      out <- 
      -( norm.const + (D/2)*log(det(t(Z)%*%Z + (sigx/siga)^2 * diag(K))) ) + 
      -1/(2*sigx^2) * tr(t(X) %*% (I-H) %*% X)
    }
    #out + dIBP(Z,a,log=T,exch=F,const=F)
    out
  } # p.x.z
  
  #p.x.z(z) # Remove this line!!!  
  p.zik.x <- function(z,i,k,log=T) { # P[z_{ik}=1|z_{-i,k}].  # exact

    zz <- z
    zz[i,k] <- 1
    a <- p.x.z(zz,log=T,with.norm.const=F) 
    zz[i,k] <- 0
    b <- p.x.z(zz,log=T,with.norm.const=F) 
    mk = sum(z[-i,k])
    p1 <- mk/N # p = Prior
    p0 <- 1 - p1

    p0 <- b+log(p0)#7Nov,2014  
    p1 <- a+log(p1)#7Nov,2014
    
    if (mk==0){
      out <- 0 # essentially doing nothing. because mk=0 => z[i,k] = 0.
    } else if (!log) { # Not Logged
      #out <- a * p / (a*p + b*(1-p))
      out <- 1 / (1+exp(p0-p1))
    } else { # Logged
      #print(paste("a:",a,"b:",b))
      #out <- a + log(p) - log(taylor.e(a)*p + taylor.e(b)*(1-p)) # HERE IS MY PROBLEM! 31 OCT
      #out <- a + log(p) - log(exp(a)*p + exp(b)*(1-p)) # HERE IS MY PROBLEM! 31 OCT
      out <- -log(1+exp(p0-p1))
    }
    
    out
  }
  
  # DOES THIS MAKE SENSE???
  # STOPPED HERE 8 NOVEMBER. START AGAIN on 10 NOVEMBER.

  # Need a way to sample alpha. Draw posterior alpha.
  #sampNewAlpha <- function(alpha=a,z,s=5) {
  #  prior <- function(x) dgamma(x,1,1,log=T)
  #  lp <- prior(seq(0,5,len=100))
  #  ll <- p.x.z(z)
  #  
  #}

  sampNewCols <- function(z,i,s=9,a=alpha[1]) { #s is the max num of columns. Avoid mh.
    prior <- function(l) dpois(l,a/N,log=T)  

    lp <- prior(0:s) # log prior
    ll <- apply(matrix(0:s),1,function(x) { 
                                col1 <- matrix(0,N,x)
                                col1[i,] <- 1
                                p.x.z(cbind(z,col1))
                                }) # log like
    lg <- lp+ll

    lpi <- apply(matrix(0:s),1,function(i) -log(sum(exp(lg-lg[i+1]))))
    Pi <- exp(lpi)
    
    #print(Pi)
    sample(0:s,1,prob=Pi)
  }

  Zs <- as.list(1:B)
  Zs[[1]] <- matrix(1,N)

  for (b in 2:B) { # B = num of iterations in Gibbs
    old.time <- Sys.time()
    z <- Zs[[b-1]]
    alpha[b] <- alpha[b-1]

    for (i in 1:N) { # iterate through all rows of Z
      K <- ncol(z)

      k <- 1
      while(k<=K) { # iterate through all columns of Z
        if (K>0) {
          #p <- p.zik.x(z,i,k,log=T) # pzik=1|x. Here
          p <- p.zik.x(z,i,k,log=F) # pzik=1|x.
          #cat(paste("\r Num of Z columns:",K))
          #Sys.sleep(.01)
          #u <- log(runif(1))
          u <- runif(1)
          if (p<=0) { # HERE!!!!!!!!!!!!!!!!!!!!!!!
            #if (p==0) print("Removed a column")
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
      #new <- rpois(1,a/N)
      new <- sampNewCols(z,i,a=alpha[b])

      if (new>0) {
        #print(paste("a/N = ", a/N ,". New Columns:",new,"ncol(z) =",ncol(z)))
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
    
    if (b %% 50 == 0) {
      gp <- b%/%50
      sink("out/Z.post.results",append=T)
        Zs[ ((gp-1)*50+1) : (gp*50) ]
      sink()  
    } else if (b==1) {
      sink("out/Z.post.results")
      sink()  
    }

    if (plotProgress) {
      n.col <- unlist(lapply(Zs[1:b],ncol))
      plot(n.col,xlab="Iteration",ylab="K+",
           main=paste0("Columns of Z ","(",b,")"),col="pink",lwd=3,type="b",pch=20)
      abline(h=mean(n.col),col="blue",lwd=3)     
      minor <- function() { plot(alpha,type="l",col="gray30") } 
      plot.in.plot(minor,"bottomright")
    }

    if (showProgress) count.down(old.time,b,B)#pb(b,B)
  } # end of gibbs

  #out <- Zs[(burn+1):B]
  out <- list("Zs"=Zs,"alpha"=alpha)
  out
}

#elapsed.time <- system.time(out <- gibbs.post(Y,a=1,B=1000,burn=0,showProgress=T,
#                                              plotProgress=T))
