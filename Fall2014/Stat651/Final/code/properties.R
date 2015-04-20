source("ibp.R")
source("rfunctions.R")

a <- 3
N <- 10
B <- 10000
Z <- lapply(as.list(1:B),function(x) rIBP(N,a))
row.sum.list <- lapply(Z,function(z) apply(z,1,sum))

row.sum <- matrix(0,0,N)
for (i in 1:B) {
  row.sum <- rbind(row.sum,unlist(row.sum.list[i]))
}

apply(row.sum,2,mean) #\bar{K1^(i)}. Property 1: a.


mean(unlist(lapply(Z,ncol))) #\bar{K+}

a * (sum(1/1:N)) # = E[K+]. Property 2.


mean(unlist(lapply(Z,sum))) #\bar{sum(Z)}
a*N # = E[sum(Z)]. Property 3.
