
samp <- function(K,d=6,N=100000) {
  
  smallSamp <- matrix(0,N,K)
  lot  <- c(rep(T,d),rep(F,100-d))
  for (i in 1:N){
    smallSamp[i,] <- sample(lot, K, replace=F)
  }
  def <- apply(smallSamp,1,sum)
  mean(def<=1)
}
theory <- ( choose(100-d,K) + choose(100-d,K-1) * choose(d,1) )/ choose(100,K)
theory
samp(51)
