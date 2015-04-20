# Pooled variance t-test (equal variance assumed)
# Returns true if the null hypothesis is rejected
t.test.pool <- function(y1,y2,alpha=0.05) {
  n1 <- length(y1)
  n2 <- length(y2)
  s1 <- sd(y1)
  s2 <- sd(y2)
  df <- n1 + n2 - 2
  s.pool <- sqrt( ( (n1-1)*s1^2 + (n2-1)*s2^2 ) / ( n1 + n2 - 2 ) )
  test.statistic <- ( mean(y1) - mean(y2) ) / ( s.pool * sqrt( 1.0/n1 + 1.0/n2 ) )
  if ( abs(test.statistic) >= qt(1-alpha/2,df) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Seperate variance t-test (Law of Large Numbers applied)
# Returns true if the null hypothesis is rejected
t.test.lln <- function(y1,y2,alpha=0.05) {
  n1 <- length(y1)
  n2 <- length(y2)
  s1 <- sd(y1)
  s2 <- sd(y2)
  #df <- n1 + n2 - 2
  test.statistic <- ( mean(y1) - mean(y2) ) / sqrt( s1^2/n1 + s2^2/n2 )
  if ( abs(test.statistic) >= qnorm(1-alpha/2) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Welch approximation t-test (degrees of freedom modified)
# Returns true if the null hypothesis is rejected
t.test.welch <- function(y1,y2,alpha=0.05) {
  n1 <- length(y1)
  n2 <- length(y2)
  s1 <- sd(y1)
  s2 <- sd(y2)
  k <- (s1^2/n1) / (s1^2/n1 + s2^2/n2)
  df <- (n1-1)*(n2-1)/((1-k)^2*(n1-1)+k^2*(n2-1))
  test.statistic <- ( mean(y1) - mean(y2) ) / sqrt( s1^2/n1 + s2^2/n2 )
  if ( abs(test.statistic) >= qt(1-alpha/2,df) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

t.test.permute <- function(y1,y2,alpha=.05) {
  n1 <- length(y1)
  n2 <- length(y2)
  test.statistic <- ( mean(y1) - mean(y2) ) / sqrt( var(y1)/n1 + var(y2)/n2 )

  test.statistic.permute <- NULL
  for (i in 1:1000) {
    y  <- sample(c(y1,y2),n1+n2)
    y1 <- y[1:n1]
    y2 <- y[-(1:n1)]
    test.statistic.permute[i] <- ( mean(y1) - mean(y2) ) / 
                              sqrt ( var(y1)/n1 + var(y2)/n2 )
  }
  quant <- quantile(test.statistic.permute,c(alpha/2,1-alpha/2))

  if ( (test.statistic <= quant[1]) || (test.statistic >= quant[2]) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }

}


# Simulation parameters
n1 <- 6
n2 <- 10
mu1 <- 0
sigma1 <- 1
mu2 <- c(0,0.5,1.0,2.0)
sigma2 <- c(0.25,1.0,2.0,5.0)
nreps <- 10000
alpha <- 0.10

# Storage for point estimates of power with row and column labels
# Note: You might be tempted to call the variable "dimnames", but that
#       would mask the function with the same name.
dnames <- list(paste("sigma2=",sigma2,sep=""),paste("mu2=",mu2,sep=""))
power.pool <- matrix(NA,nrow=length(sigma2),ncol=length(mu2),dimnames=dnames)
power.lln <- matrix(NA,nrow=length(sigma2),ncol=length(mu2),dimnames=dnames)
power.welch <- matrix(NA,nrow=length(sigma2),ncol=length(mu2),dimnames=dnames)
power.permute <- matrix(NA,nrow=length(sigma2),ncol=length(mu2),dimnames=dnames)

# Storage for confidence intervals of power
dnames <- list(paste("sigma2=",sigma2,sep=""),rep(paste("mu2=",mu2,sep=""),each=2))
ci.pool <- matrix(NA,nrow=length(sigma2),ncol=2*length(mu2),dimnames=dnames)
ci.lln <- matrix(NA,nrow=length(sigma2),ncol=2*length(mu2),dimnames=dnames)
ci.welch <- matrix(NA,nrow=length(sigma2),ncol=2*length(mu2),dimnames=dnames)
ci.permute <- matrix(NA,nrow=length(sigma2),ncol=2*length(mu2),dimnames=dnames)

library(foreach)
library(doMC)
registerDoMC(32)
useParallel <- T


for ( i in 1:length(sigma2) ) {
  for ( j in 1:length(mu2) ) {
    results <- matrix(NA,nrow=nreps,ncol=4)

    engine <- function() {
      y1 <- rnorm(n1,mean=mu1,sd=sigma1)
      y2 <- rnorm(n2,mean=mu2[j],sd=sigma2[i])
      c(t.test.pool(y1,y2,alpha),t.test.lln(y1,y2,alpha),
        t.test.welch(y1,y2,alpha),t.test.permute(y1,y2,alpha))
    }
    if (!useParallel) {
      for (k in 1:nreps){
        results[k,] <- engine()
      }
    } else {
      results <- foreach(k=1:nreps,.combine=rbind) %dopar% engine()
    }

    # Point estimates
    power.pool[i,j]  <- mean(results[,1])
    power.lln[i,j]   <- mean(results[,2])
    power.welch[i,j] <- mean(results[,3])
    power.permute[i,j] <- mean(results[,4])
    # 95% confidence interval
    ci.pool[i,rep(2*(j-1),2)+c(1,2)]  <- power.pool[i,j]  + c(-1,1) * qnorm(1-0.05/2) * sqrt( var(results[,1])/nreps )
    ci.lln[i,rep(2*(j-1),2)+c(1,2)]   <- power.lln[i,j]   + c(-1,1) * qnorm(1-0.05/2) * sqrt( var(results[,2])/nreps )
    ci.welch[i,rep(2*(j-1),2)+c(1,2)] <- power.welch[i,j] + c(-1,1) * qnorm(1-0.05/2) * sqrt( var(results[,3])/nreps )
    ci.permute[i,rep(2*(j-1),2)+c(1,2)]  <- power.permute[i,j]  + c(-1,1) * qnorm(1-0.05/2) * sqrt( var(results[,4])/nreps )

    print(c(i,j))
  }
}

power.pool
power.lln
power.welch
power.permute

ci.pool
ci.lln
ci.welch
ci.permute
##############################################################################


