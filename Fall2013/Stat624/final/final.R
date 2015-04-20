logit <- function(x){
  log(x/(1-x))
}

f <- function(x,mu,sig2){
  ans <- ifelse((x>0) & (x<1),1/(sqrt(sig2*2*pi)) * exp(-(logit(x)-mu)^2/(2*sig2)) / (x*(1-x)),0)
}

mean.via.direct.sampling <- function(mu,sigma2,B){
  y <- rnorm(B,mu,sqrt(sigma2))
  x <- exp(y) / (1+exp(y))
  CI <- qnorm(c(.025,.975),mean(x),sd(x)/sqrt(B))
  result <- matrix(c(mean(x),CI[1],CI[2]),nrow=1)
  colnames(result) <- c("Estimate","CI Lower","CI Upper")
  result
}

mean.via.rejection.sampling <- function(mu,sigma2,B){
  
  rejection.sampler <- function(log.density.target,log.density.envelope,
                                sample.envelope,alpha,sample.size) {

    n <- sample.size
    x <- matrix(0,n,1)
    if ( is.vector(sample.envelope()) ) {
      x <- matrix(0,n,length(sample.envelope()))
    }

    my.sampler <- function(x){
      out <- sample.envelope()
      while ( log.density.target(out)-(log.density.envelope(out)-log(alpha)) <= log(runif(1)) ){
        out <- sample.envelope()
      }
      out
    }

    return (t(apply(x,1,my.sampler)))
  }

  #mu <- 2; sigma2 <- 1
  #curve(f(x,mu,sigma2),from=0,to=1)
  h <- function(x) -logit(x) + sigma2*(2*x-1) + mu
  mode <- uniroot(h,c(0,1))$root
  #abline(v=mode$root)

  log.f <- function(x){
    -(logit(x)-mu)^2/(2*sigma2) - log(x*(1-x)) 
  }
  
  log.g <- function(x){
    0
  }

  sampler <- function() runif(1)
  alpha <- exp(log.g()-log.f(mode))
  x <- rejection.sampler(log.f,log.g,sampler,alpha,B)

  #CI <- mean(x) + c(-1,1)*qnorm(.975)*sd(x)/sqrt(B)
  CI <- qnorm(c(.025,.975),mean(x),sd(x)/sqrt(B))
  result <- matrix(c(mean(x),CI[1],CI[2]),nrow=1)
  colnames(result) <- c("Estimate","CI Lower","CI Upper")
  result
}

mean.via.importance.sampling <- function(mu,sigma2,B){
  h <- function(x) x*f(x,mu,sigma2)
  g <- function(x) 1 # because the pdf for U(0,1) = 1
  u <- runif(B)
  x <- h(u)/g(u)
  #CI <- mean(x) + c(-1,1)*qnorm(.975)*sd(x)/sqrt(B)
  CI <- qnorm(c(.025,.975),mean(x),sd(x)/sqrt(B))
  result <- matrix(c(mean(x),CI[1],CI[2]),nrow=1)
  colnames(result) <- c("Estimate","CI Lower","CI Upper")
  result
}

mean.via.mcmc <- function(mu,sigma2,B){
  out <- NULL
  cs <- .35
  cnt <- 0
  out[1] <- .5
  
  lf <- function(x) -(logit(x)-mu)^2/(2*sigma2) - log(x*(1-x)) 

  for (i in 2:B){
    out[i] <- out[i-1]
    cand <- rnorm(1,out[i],cs)
    if( (cand>0) & (cand<1) ){
      a <- lf(cand) - lf(out[i])
      if (a > log(runif(1))){
        out[i] <- cand
        cnt <- cnt + 1
      }  
    }
  }

  #plot(out[-c(1:round(B*.2))],type='l')
  #print(cnt/B)
  
  x <- out[-c(1:round(B*.2))]
  #CI <- mean(x) + c(-1,1)*qnorm(.975)*sd(x)/sqrt(length(x))
  CI <- qnorm(c(.025,.975),mean(x),sd(x)/sqrt(B))
  result <- matrix(c(mean(x),CI[1],CI[2]),nrow=1)
  colnames(result) <- c("Estimate","CI Lower","CI Upper")
  result
}

mean.via.simpsons.rule <- function(mu,sigma2,M){
  fun <- function(x) x*f(x,mu,sigma2)
  a <- 0; b <- 1
  h <- (b-a)/(2*M)
  x.v <- seq(a,b,by=h)
  f.v <- sapply(x.v,fun)
  h * (f.v[1]/2+sum(f.v[2:(2*M)])+f.v[2*M+1]/2)
}

