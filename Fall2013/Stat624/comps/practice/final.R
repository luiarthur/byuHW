L <- 1000000

logit <- function(x){
  log(x/(1-x))
}

f <- function(x,mu,sig2){
  ans <- ifelse((x>0) & (x<1),1/(sqrt(sig2*2*pi)) * exp(-(logit(x)-mu)^2/(2*sig2)) / (x*(1-x)),0)
  ans
}

display <- function(x,b) {
  CI <- matrix(c(mean(x),qnorm(c(.025,.975),mean(x),sd(x)/sqrt(b))),1)
  colnames(CI) <- c("Estimate","CI.Lower","CI.Upper")
  CI
}

mean.via.direct.sampling <- function(mu,sigma2,B){
  y <- rnorm(B,mu,sqrt(sigma2))
  x <- exp(y) / (1+exp(y))
  display(x,B)
}
#mean.via.direct.sampling(.5,2,L)


mean.via.importance.sampling <- function(mu,sigma2,B){
  x <- runif(B) # X ~ U[0,1]
  h <- x * f(x,mu,sigma2)
  g <- 1 # the pdf for U[0,1]
  display(h/g,B)
}
#mean.via.importance.sampling(.5,2,L)

mean.via.mcmc <- function(mu,sigma2,B){
  out <- NULL
  out[1] <- .5
  cs <- .73
  acc <- 0
  burn <- round(B*.9)

  lf <- function(x) {
    -(logit(x)-mu)^2/(2*sigma2) - log(x) - log(1-x)
  }
  
  pb <- txtProgressBar(2,B,width=30,style=3)

  for (i in 2:B) {
    out[i] <- out[i-1]
    cand <- rnorm(1,out[i],cs)
    
    if ( (cand>0) & (cand<1) ) {
      z <- lf(cand) - lf(out[i])
      if (z > log(runif(1))) {
        out[i] <- cand
        if (i > burn) acc <- acc + 1
      }
    }
    
    setTxtProgressBar(pb,i)
  }
  
  out <- tail(out,burn)
  #print(acc/(B-burn))
  #plot(out)
  #plot(density(out))
  display(out,B-burn)
}
#temp <- mean.via.mcmc(.5,2,L)


mean.via.rejection.sampling <- function(mu,sigma2,B){
  # Uniform Envelop
  U <- runif(B)
  X <- runif(B)
  mode.fn <- function(x) sigma2*(2*x-1)+mu -logit(x)
  mode <- uniroot(mode.fn,c(0,1))$root 
  mx <-  f(mode,mu,sigma2)
  sim <- f(X,mu,sigma2) / mx
  draws <- X[which(sim > U)]
  display(draws,length(draws))
}

#mean.via.rejection.sampling(.5,2,L)
