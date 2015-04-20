test <- T

# Prior Parameteres:
pr.mu <- 0
pr.s  <- 3

dat <- rnorm(200,5,3)

log.prior <- function(x) {
  dnorm(x,pr.mu,pr.s,log=T) 
}

llik <- function(mu) {
  ifelse(test,0,sum(dnorm(dat,mu,3,log=T)))
}

llik.pr <- function(mu) {
  llik(mu) + log.prior(mu)
}

mh <- function(cs=3,B=10000,init=pr.mu,burn=B*.1) {
  pb <- txtProgressBar(min=2,max=B,style=3,width=30)
  out <- NULL
  out[1] <- init
  cnt <- 0

  for (i in 2:B) {
    setTxtProgressBar(pb,i)
    out[i] <- out[i-1]
    cand <- rnorm(1,out[i],cs)

    r <- llik.pr(cand) - llik.pr(out[i])
    if (r > log(runif(1))) {
      out[i] <- cand
      if (i > burn) cnt <- cnt + 1
    }
  }

  close(pb)
  print(cnt/(B-burn))
  tail(out,B-burn)
}

plot.diag <- function(x) {
  par(mfrow=c(2,1))
    plot(x,type="l")
    plot(density(x),col="blue")
    curve(dnorm(x,pr.mu,pr.s),from=min(out),to=max(out),col="red",add=T)
    legend("topleft",legend=c("Prior","Posterior"),col=c("red","blue"),lwd=3)
  par(mfrow=c(1,1))
}

boot <- function(data,fn,N=10000) {
  n <- length(data)
  one.it <- function(x) {
    d <- sample(data,n,replace=T)
    fn(d)
  }
  res <- apply(matrix(1:N),1,one.it)
  est.CI <- matrix(c(mean(res),qnorm(c(.025,.975),mean(res),sd(res))),1)
  colnames(est.CI) <- c("Estimate","CI.Lower","CI.Upper")
  est.CI
}

out <- mh(cs=10)
plot.diag(out)

result <- t(sapply(list(mean,sd),function(x) boot(out,x)))
result <- cbind(result,matrix(c(pr.mu,pr.s)))
contained <- apply(result,1,function(x) ifelse(x[2] < x[4] & x[4] < x[3],1,0))
result <- cbind(result,contained)

colnames(result) <- c("Estimate","CI.Lower","CI.Upper","True","Contained")
rownames(result) <- c("Mean","sd")


result

# So, if you set the likelihood to 0, you are sampling from the prior.
