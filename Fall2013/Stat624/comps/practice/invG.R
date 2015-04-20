a <- 3
b <- 4
B <- 10^5

f <- function(x) {
  b^a / gamma(a) * x^(-a-1) * exp(-b/x)
}

display <- function(x,L=B) {
  CI <- matrix(c(mean(x),qnorm(c(.025,.975),mean(x),sd(x)/sqrt(L))),nrow=1)
  colnames(CI) <- c("Estimates","CI Lower","CI Upper")
  CI
}


dr <- function() {
  Y <- 1/rgamma(B,a,scale=1/b)
  display(Y)
}
#dr()

im <- function() {
  inf <- 10^3
  X <- runif(B,0,inf)
  g <- 1/inf #function(x) dgamma(x,2,2)
  h <- X * f(X)
  X <- h/g
  display(X)
}
#im()

B <- 100000
mc <- function() {
  out <- NULL
  cs <- 2
  acc <- 0
  burn <- round(.9*B)
  out[1] <- 1

  pb <- txtProgressBar(2,B,style=3)
  
  for (i in 2:B) {
    out[i] <- out[i-1]
    cand <- rnorm(1,out[i],cs)
    
    if(cand>0){
      r <- f(cand) / f(out[i])
      if(r > runif(1)) {
        out[i] <- cand
        if (i > burn) acc <- acc + 1
      }
    }
    
    setTxtProgressBar(pb,i)
  }
  #close(pb)

  out <- tail(out,B-burn)
  #print(acc/length(out))
  #plot(out)
  display(out)
}
#mc()

rj <- function() {
  mode <- b/(a+1)
  inf <- 10^3
  X <- runif(B,0,inf)
  U <- runif(B)
  mx <- f(mode)
  sims <- f(X) / mx 
  draws <- X[sims > U]
  display(draws)
}
#rj()


library(foreach)
library(doMC)
registerDoMC(4)

fn <- list(dr,im,mc,rj)
M <- foreach(i=1:4,.combine=rbind) %dopar% fn[[i]]()
rownames(M) <- c("Direct","Importance","MCMC","Rejection")

M
