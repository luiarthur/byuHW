n <- 3
results <- NULL
N <- 10000
for (i in 1:N){
  y <- rgamma(n,3,5)
  results[i] <- (max(y) > qgamma(.5,3,5))
}
p <- mean(results)
CI <- p + c(-1,1) * qnorm(.975) * sqrt(p * (1-p) / N) 
ans <- matrix(c(p, CI, 1-.5^n),nrow=1)
colnames(ans) <- c('Empirical','CI.low','CI.Hi','Theoretical')

ans
