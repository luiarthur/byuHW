mann.whitney <- function(sample1,sample2,alternative=c("greater","less")[1]){

  Np <- 5000
  U <- NULL
  alpha <- .05

  n1 <- length(sample1)
  n2 <- length(sample2)

  s1 <- cbind(sample1,rep(1,n1))
  s2 <- cbind(sample2,rep(2,n2))

  s <- rbind(s1,s2)
  so <- s[order(s[,1]),]

  R1 <- sum(which(so[,2] ==1))
  U1 <- R1 - n1*(n1+1)/2 

  for (i in 1:Np){

      #browser()
      s1 <- rnorm(n1)
      s2 <- rnorm(n2)

      s1 <- cbind(s1,rep(1,n1))
      s2 <- cbind(s2,rep(2,n2))

      s <- rbind(s1,s2)
      so <- s[order(s[,1]),]

      R1 <- sum(which(so[,2] ==1))
      U[i] <- R1 - n1*(n1+1)/2 

  }

  mc.pVal <- ifelse(alternative == 'greater', mean(U>=U1), mean(U<U1) )

  CI <- mean(mc.pVal) -c(-1,1) * qnorm(alpha/2)* sqrt(mc.pVal*(1-mc.pVal)/Np)
  c(U1, mc.pVal, CI[1], CI[2])
}


# To commit: (1) Ctrl-V, (2) Shift i, (3) Type '#', (4) ESC, (5) Wait a second
#

samp1 <- rnorm(20) 
samp2 <- rnorm(20)

verify <- function(){

  mw <- mann.whitney(samp1, samp2, "greater")

  cat("The Test Statistic U is:\t", mw[1], '\n')
  cat("The Monte Carlo p-value is:\t", mw[2], "(",mw[3],",",mw[4],")", '\n')
}


verify()
wilcox.test(samp1,samp2,'greater')
