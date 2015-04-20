# Transforming parameters (which are not random) is fine.
# Transforming R.V.'s can be problematic.
# To interpret beta, consider a plot of Y ~ Xj, HOLDING ALL OTHER X's CONSTANT!

library(foreach)
library(doMC)
registerDoMC(16)

# BIC function:
bic <- function(y,x,b){
  b <- as.matrix(b)
  p <- length(b) - 1 
  n <- length(y)
  e <- y - x%*%b
  rss <- t(e)%*%e 
  -2*log(rss) + p*log(n)
}

bayes.fwd <- function(y=crash$Fatal,x=model.matrix(y~.,data=crash),p=ncol(x)-1){
  M <- list(); length(M) <- p
  Bic <- NULL
  ind <- 1 
  pool <- 2:(p+1)
  for (k in 1:p){
    BIC <- NULL

    # Sequential
    #for (i in pool){
    #  res <- bayes.probit(100,y,x[,c(ind,i)]); cat(paste(k,i,collapse=""))
    #  b <- apply(res,2,mean)
    #  BIC[k] <- bic(y,x[,c(ind,i)],b)
    #}

    # Parallel
    one.it <- function(i) {  
      res <- bayes.probit(100,y,x[,c(ind,i)]); cat(paste(k,i,collapse=""))
      b <- apply(res,2,mean)
      BIC[k] <- bic(y,x[,c(ind,i)],b)
    }
    BIC <- as.vector(foreach(i=pool,.combine=rbind) %dopar% one.it(i))

    new <- pool[which.min(BIC)]
    ind <- c(ind,new)
    Bic[k] <- BIC[new]
    #M[[k+1]] <- ind
    M[[k]] <- ind
    pool <- pool[-which.min(BIC)]
  }
  list("bic"=Bic,"M"=M,"best.set"=M[which.min(Bic)])
}

#mod.fwd <- bayes.fwd()

# Libraries Need:
options("scipen"=8)
options("width"=120)
library(truncnorm) # for rtruncnorm
pc <- function(X,M) t(eigen(cov(X))$vectors[,1:M])
mvrnorm <- function(M,S,n=nrow(S))  M + t(chol(S)) %*% rnorm(n)

# Data Cleaning:
crash <- read.csv("../Data/crash.csv")
crash <- crash[which(crash$Hour<=23),] # remove the 99th Hour. 23 of them.
crash <- crash[-which(crash$Mod_year>2013),] # remove the 9999 model year. 4 of them. Come back to this.
crash$Mod_year[which(crash$Mod_year<1987)] <- 1986 # Since there are only a few observations less than 1987, I grouped them together.
crash <- crash[,-2] # There is only one year, so remove the Year column
#colnames(crash)
#str(crash)

# Probit Bayesian Analysis:
# Yi = I{Zi > 0} ~ Bern(pi)
# pi = F(Xb), where F is the CDF for Normal(0,1)
# Zi ~ N(Xb,1), the 1 is arbitrary and works for this bayesian analysis
# P(Yi=1) = P(Zi > 0) = F(Xb)

updateZ <- function(x,y,b){
  z <- numeric(length(y))
  z[y==1] <- rtruncnorm(sum(y==1),a=0,b=Inf,mean=x[y==1,] %*% b,sd=1)
  z[y==0] <- rtruncnorm(sum(y==0),a=-Inf,b=0,mean=x[y==0,] %*% b,sd=1)
  z
}


bayes.probit <- function(B=10000,Y=crash$Fatal,X=model.matrix(Y~.,data=crash),...){
  
  XtXi <- solve(t(X)%*%X)
  Xt <- t(X)

  # Initialize Parameters
  z <- 1
  beta <- matrix(0,B,ncol(X))
  #######################

  pb <- txtProgressBar(min=2,max=B,width=30,...) # Define Progress Bar
  for (i in 2:B){
    #Updates:
    setTxtProgressBar(pb, i) # Update Progress Bar
    z <- updateZ(X,Y,beta[i,])
    beta[i,] <- mvrnorm(XtXi %*% Xt%*%z, XtXi) 
  }
  close(pb) # Close Progress Bar

  beta[-c(1:B%/%10),]
}

one.sim <- function(BB=10000,cv=F,size=0,PCR=F,M=0,thresh=.5) {

  # Set Y and X
  testI <- sample(1:nrow(crash),size,repl=F)
  if (!cv) testI <- -(1:nrow(crash))

  # Training Set
  Y <- crash$Fatal[-testI]
  X <- model.matrix(Fatal ~ .,data=crash[-testI,])
  n <- nrow(X)
  
  PSI <- NULL
  # Testing PCR:#######
  if (PCR) {
    PSI <- pc(scale(X[,-1]),M)
    Z <- X[,-1] %*% t(PSI)
    Z <- cbind(1,Z)
    X <- Z
  }
  #####################

  comp.time <- system.time(result <- bayes.probit(B=BB,Y,X,style=3))
  mean.beta <- apply(result,2,mean)
  se.beta <- apply(result,2,sd)
  s2 <- t(Y-X%*%mean.beta) %*% (Y-X%*%mean.beta) / (n-length(mean.beta))
  if (PCR) {
    mean.beta <- rbind(mean.beta[1],t(PSI) %*% mean.beta[-1])
    se.beta <- c(se.beta[1],
                 sqrt(s2 * diag(t(PSI) %*% solve(t(X[,-1])%*%X[,-1]) %*% PSI)))
  }

  MS <- cbind(mean.beta,se.beta)
  get.ci <- function(ms) t(qnorm(c(.025,.975),ms[1],ms[2]))
  if (PCR) get.ci <- function(ms) t(ms[1]+qt(c(.025,.975),n-M-1)*ms[2])
  MS.CI <- as.data.frame(cbind(MS,t(apply(cbind(MS),1,get.ci))))
  signif <- ifelse(!MS.CI[,3] <= 0 & 0 <= MS.CI[,4],"*","")
  MS.CI <- cbind(MS.CI,signif)
  colnames(MS.CI) <- c("Estimate","Std.Err","CI.Lower","CI.Upper","Sig")
  #rownames(MS.CI) <- paste("beta.",0:(ncol(X)-1),sep="")
  rownames(MS.CI) <- colnames(X)

  signif.beta <- MS.CI[which(MS.CI$Sig=="*"),]
  # Now need to do PCR

  if (!cv) { 
    BIC <- bic(Y,X,MS.CI$Est)
    list("MS.CI"=MS.CI,"signif"=signif.beta,"X"=X,"Y"=Y,"s2"=s2,"bic"=BIC,"beta"=result)
  } else {

    # Low Accuracy: 25% Error Rate
    sig <- which(signif=="*")
    xb <- model.matrix(Fatal ~ ., data=crash[testI,])[,sig] %*% mean.beta[sig]
    pred <- ifelse(pnorm(xb)>thresh,1,0)

    # High Accuracy: 17% Error Rate
    #xb <- model.matrix(Fatal ~ ., data= crash[testI,]) %*% mean.beta # High Acc.
    #pred <- ifelse(xb>thresh,1,0)

    true <- crash$Fatal[testI]
    typy <- sum(true==1 & pred==1)
    tnpn <- sum(true==0 & pred==0)
    sens <- typy / sum(true==1)
    spec <- tnpn / sum(true==0)
    err.rate <- mean(true!=pred)

    list("MS.CI"=MS.CI,"signif"=signif.beta,"X"=X,"Y"=Y,"Sens"=sens,"Spec"=spec,
         "err.rate"=err.rate,"s2"=s2,"beta"=result)
  }
}  

#temp.1 <- one.sim(1000)

# Get Sens and Spec:
N <- 100
thresh <- 1:N/N
f <- function(i) {print(i); one.sim(100,cv=T,size=100,thresh=thresh[i])}
result <- foreach(j=1:N,.errorhandling="remove") %dopar% f(j)
sens <- sapply(result,function(x) x$Sens)
spec <- sapply(result,function(x) x$Spec)
err  <- sapply(result,function(x) x$err.rate)
write.table(cbind(sens,spec),"out/results.txt",quote=F,row=F)

plot.trace <- function(out,o=4) {
  par(mfrow=c(o,1))    
    for (i in 1:63){
      #plot(temp.1$beta[,i],type='l',main=paste(expression(beta),i-1))
      plot(out[,i],type='l',main=paste(expression(beta),i-1))
      if (i%%o==0) {locator(1)}
    }
  par(mfrow=c(1,1))
}

plot.trace(result[[1]]$beta,4)


# Plot ROC
plot(1-spec,sens,xlim=c(0,1),ylim=c(0,1),col="blue",cex=.5,main="ROC"); abline(0,1)
#lines(lowess(1-spec,sens))
mod <- lm(sens ~ I(log(1.0001-spec)))
h <- function(x,beta) beta[1] + beta[2] * log(x)
curve(h(x,beta=mod$coef),fr=.0001,to=.97,add=T,col='blue',lwd=2)
mean(err)


# Find Best Threshold:
opt.thresh <- which.min(sapply(result,function(x) x$err.rate)) / 100



########################################################################
# Get M for PCR:
#m <- 1:60
#g <- function(i) {print(i); one.sim(1000,cv=T,size=100,M=m[i],PCR=T)}
#result <- foreach(i=m,.errorhandling="remove") %dopar% g(i)
#sens <- sapply(result,function(x) x$Sens)
#spec <- sapply(result,function(x) x$Spec)
#write.table(cbind(sens,spec),"out/results.txt",quote=F,row=F)
#plot(1-spec,sens,xlim=c(0,1),ylim=c(0,1),col="red",cex=1,pch=20); abline(0,1)
#plot(identify(1-spec,sens),xlim=c(0,1),ylim=c(0,1),col="red",cex=1,pch=20); abline(0,1)
## Best to have M around 45

# Compare PCR and Normal
#result.pcr <- one.sim(1000,M=5,PCR=T)
#result.nor <- one.sim(1000)
#cbind(result.pcr$MS.CI$Est,result.nor$MS.CI$Est)
#result.pcr$signif
#result.nor$signif

