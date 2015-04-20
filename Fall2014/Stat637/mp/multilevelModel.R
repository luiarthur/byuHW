# Paper:
#http://www.stat.columbia.edu/~gelman/research/published/multi2.pdf
# Data taken from:
#http://www.unc.edu/courses/2007spring/enst/562/001/docs/assignments/assign10.htm

det <- function(x,log=F) {
  out <- 0
  if (!log) {
    out <- det(x)
  } else {
    out <- unlist(determinant(x,log=T))[1]
  }
  out
}

radon <- read.csv("radon.txt",header=T)
city <- read.csv("cty.txt",header=T)

rad <- radon[which(radon$state=="MN"),c("activity","basement","county")]
rad <- rad[which(sapply(rad$basement,as.character) > ""),]
rad <- rad[which(rad$activity > 0),]
ctyInfo <- city[which(city$st=="MN"),c("cty","Uppm")]
dat <- merge(rad,ctyInfo,by.x="county",by.y="cty")

n <- nrow(dat)
k <- length(unique(dat$county))
y <- log(dat$activity)
x <- ifelse(dat$basement=="Y",1,0)
counties <- unique(dat$county)
G <- matrix(0,n,k)
u <- unique(dat$Uppm)
for (i in 1:n) G[i,which(counties==dat$county[i])] <- 1 
Gu <- G%*%u
G2 <- G%*%t(G)

ll <- function(g0,g1,b,sy2,sa2) {# Proportional To. Log likelihood.
  m <- g0 + g1*Gu + b*x
  S <- diag(sy2,n) + sa2*G2
  out <- -.5*(t(y-m) %*% solve(S) %*% (y-m) - det(S,log=T))
  out
}

#Priors:
# a ~ N_k(0,100I)
# b ~ N(0,100)
# g_0 ~ U(0,100)
# g_1 ~ U(0,100)
# sy2 ~ U(0,100)
# sa2 ~ U(0,100)

mh <- function(B=1e4,candSigA=diag(100,k)) {
  A <- matrix(0,B,k)
  bb <- rep(0,B)
  G0 <- rep(1,B)
  G1 <- rep(1,B)
  Sy2 <- rep(1,B)
  Sa2 <- rep(1,B)
  
  for (i in 2:B) {
    A[i,] <- A[i-1,]
    bb[i] <- bb[i-1]
    G0[i] <- G0[i-1]
    G1[i] <- G1[i-1]
    Sy2[i] <- Sy2[i-1]
    Sa2[i] <- Sa2[i-1]

    # Update A:
    cand <- mvrnorm(k,A[i,],candSigA)
    r <- ll() + lp 
  }

}
