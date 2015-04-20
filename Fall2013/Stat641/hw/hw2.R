N <- 10^6
n <- 5

simulate <- function(x){
  balls <- rmultinom(N, x, rep(1/x,x))
  balls 
}

countExactOneEmpty <- function(sim){
  count <- 0
  for (i in 1:(dim(sim)[2])){
    if ( length(which(sim[,i] == 0)) == 1){
      count = count + 1
    }
  }
  count
}
#sims <- simulate(n)
#countExactOneEmpty(sims) / N


cat('n', '\t', 'Empirical', '\t', 'Theoretical', '\n')
for ( i in 5:10 ) {
  cat(i, '\t', countExactOneEmpty(simulate(i))/N,
         '\t', factorial(i) * choose(i,2) / i^i,'\n')
}

