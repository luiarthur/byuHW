updateMu <- function(x,v,mu,a){
  n <- length(x)
  return( c((mu*v+sum(x)*a)/(v+n*a),(a*v)/(v+n*a)) )
}
updateSS <- function(x,mu,sh,sc){
  n <- length(x)
  sx <- sum((x-mu)^2)
  return (c(n/2+sh, 2*sc/(sx*sc+2)))
}

gibbsN <- function(x,mumu,ssmu,a,b,N=100){
 
  out <- matrix(0,N,2)
  out[1,1] <- mean(x)
  out[1,2] <- var(x)
 
  for (j in 2:N){
  postmu <- updateMu(x,out[j-1,2],mumu,ssmu)
  out[j,1] <- rnorm(1,postmu[1],sqrt(postmu[2]))
  postSS <- updateSS(x,out[j,1],a,b)
  out[j,2] <- 1/rgamma(1,postSS[1],scale=postSS[2])
  }

  return(out)
}

x <- rnorm(1000,5,4)
out <- gibbsN(x,3,8,2,5,50000)
apply(out,2,mean)

# This works!
# Note that the mean and VARIANCE (not sd) are returned.
