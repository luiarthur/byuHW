rm(list=ls())
source('frechet.R')

pow <- function(a,b){
  a^b
}

#temp1 <- rfrechet(1000000,6,5,3)
#theoretical.stat.frechet(6,5,3)
#mean(temp1); var(temp1)
#
#temp2 <- rfrechet(1000000,6) * 3 + 5
#mean(temp2); var(temp2)

# test c code:

lL <- function(x, n, v){
  a <- v[1]; m <- v[2]; s <- v[3]
  z <- (x-m)/s
  n*log(a) - n*log(s) - sum((1+a)*log(z)+pow(z,-a))
}
lLa <- function (x, n, a, m, s){
  z <- (x-m)/s
  n/a + sum(log(z) * (pow(z,-a) - 1))
}
lLm <- function(x, n, a, m, s){
  z <- (x-m)/s
  sum((a + 1 - a*pow(z,-a)) / (z*s))
}
lLs <- function (x, n, a, m, s){
  z <- (x-m)/s
  -n/s + sum(a+1 - a*pow(z,-a))/s
}

# double derivatives: 6 unique ones 
lfaa <- function (x, n, a, m, s){
  z <- (x-m)/s
  -n/(a^2) - sum(log(z)^2 * pow(z,-a));
}

lfmm <- function (x, n, a, m, s){
  z <- (x-m)/s
  sum( (a+1)/pow(z*s,2) * (1-a*pow(z,-a)) )
}

lfss <- function (x, n, a, m, s){
  z <- (x-m)/s
  pow(s,-2) * ( n + sum(-(a+1) - a*(a-1)*pow(z,-a)) )
}

lfam <- function (x, n, a, m, s){
  z <- (x-m)/s
  sum(((a * log(z) - 1) * pow(z,-a) + 1) / (z*s))
  #sum(((a * log(z) - 1 + pow(z,a)) * pow(z,-a)) / (z*s))
}

lfms <- function (x, n, a, m, s){
  z <- (x-m)/s
  -a^2 * pow(s,a-1) * sum(pow(z*s,-a-1)); 
}

lfas <- function (x, n, a, m, s){
  z <- (x-m)/s
  sum( (a*log(z) + pow(z,a) - 1) / (s*pow(z,a)) )
}

hessian <- function(x,n,a,m,s){
  aa <- lfaa(x,n,a,m,s)
  am <- lfam(x,n,a,m,s)
  as <- lfas(x,n,a,m,s)
  ms <- lfms(x,n,a,m,s)
  mm <- lfmm(x,n,a,m,s)
  ss <- lfss(x,n,a,m,s)

  matrix(c(aa,am,as,
           am,mm,ms,
           as,ms,ss),3,3,byrow=T)
}

gradient <- function(x,n,a,m,s){
  matrix(c(lLa(x,n,a,m,s),
           lLm(x,n,a,m,s),
           lLs(x,n,a,m,s)),ncol=1)
}

newton <- function(x,n,v,r=1,max.it=200,tol=.1){
  g <- gradient(x,n,v[1],v[2],v[3])
  H <- hessian(x,n,v[1],v[2],v[3])
  print(H)
  N <- 1
  while ((N<=max.it) & (sqrt(sum(g^2))>tol) ){ 
    #print(paste("Iteration: ",N))
    v <- v - solve(H) %*%  g
    #v <- v - r*solve(H,g)
    a <- v[1]; m <- v[2]; s <- v[3]
    H <- hessian(x,n,a,m,s)
    g <- gradient(x,n,a,m,s)
    N <- N + 1  
    #print(H)
  }
  print(t(v))
}

n <- 100000; a <- 8; m <- 4; s <- 1
x <- rfrechet(n,a,m,s)
v <- c(a,m,s)

dyn.load("c/frechet.so")
.C("newton",as.double(x),as.integer(n),estimates=as.double(v))$estimates

my.lL <- function(p) -lL(x,n,p)
newton(x,n,v) 


# This tells me that at least my loglikelihood is correct
newton(x,n,v) 
est <- nlm(my.lL, init)$estimate
lines(density(rfrechet(10000,est[1],est[2],est[3])), col='blue',lwd=3)

lL2 <- function(p){
  a <- p[1]; m <- p[2]; s <- p[3]
  -sum(log(dfrechet(x,a,m,s)))
}
est <- nlm(lL2,c(4,3,3))$estimate
lines(density(rfrechet(10000,est[1],est[2],est[3])),lwd=3,col='blue')


# Metropolis:
dyn.load("c/frechet.so")
n <- 100;  N <- 100000
#a <- 5; m <- 11; s <- 6
#x <- rfrechet(n,a,m,s)
#init <- c(a,m,s)

#x <- carbone
#n <- length(x)
#init <- c(39,-33,38)

data(quakes)
x <- quakes$mag
init <- c(19,-1,5)
min <- min(x)
outA <- rep(0,N)
outM <- rep(0,N)
outS <- rep(0,N)
#par(mfrow=c(1,1))
cs <- c(3,.10,.10)
Cdata <- .C("mig",as.double(x),as.integer(n),as.integer(N),as.double(min),
                  A=as.double(outA), M=as.double(outM), S=as.double(outS),
                  as.double(init),as.double(cs))
A <- Cdata$A; plot(A,type='l'); mA <- mean(A[90000:100000])
M <- Cdata$M; plot(M,type='l'); mM <- mean(M[90000:100000])
S <- Cdata$S; plot(S,type='l'); mS <- mean(S[90000:100000])

#plot(density(x),lwd=3,ylim=c(0,1.2))
hist(x,prob=T,ylim=c(0,1.2))
curve(dfrechet(x,mA,mM,mS),from=4,to=6.5,add=T,col='red',lwd=3)

#mh <- function(x,init,min=min(x),cs=c(1,.3,.1),N=1000){
#  
#  out <- matrix(0,N,3)
#  out[1,] <- init 
#  count <- c(0,0,0)
#  n <- length(x)
#
#  for (i in 2:N){
#    out[i,] <- out[i-1,]
#
#    lp <- function(t,k,v){
#      v[k] <- t
#      a <- v[1]; m <- v[2]; s <- v[3]
#      z <- (x-m)/s
#      if (k==1){
#        n*log(a) -  sum( (a) * log(z) + z^-a) 
#        #(mean(x)-1)*log(a) - a/1 #Prior: logdgamma centered at mean(x)
#        #(5-1)*log(5) - 5/1 #Prior: logdgamma centered at mean(x)
#      } else if (k==2) {  
#        - sum((1+a)*log(z) + z^-a)  
#        #-((m-min(x)*.8)/1)^2 / 2 #Prior: logdnorm centered at min(x)
#        #-((m-0)/1)^2 / 2 #Prior: logdnorm centered at mean(x)
#      } else {  
#        -n*log(s) -sum( (1+a)*log(z) + z^-a ) 
#        #(sd(x)-1)*log(s) - s/1 #Prior: logdgamma centered at sd(x)
#        #(1-1)*log(s) - 1/1 #Prior: logdgamma centered at sd(x)
#      }  
#    }
#
#    for (j in 1:3){
#      cand <- rnorm(1,out[i,j],cs[j])
#      good <- ifelse(j==2,cand<min, cand>0)
#      if (good){
#        v <- out[i,]
#        r <- lp(cand,j,v) - lp(out[i,j],j,v)
#        if (r > log(runif(1))){
#          out[i,j] <- cand
#          count[j] <- count[j] + 1
#        }
#      }
#    }
#  } #End of MH Loop
#
#  print(count/N)
#  out
#}
#
#a <- 17; m <- 5; s <- 2
#x <- rfrechet(1000,a,m,s)
#v <- c(mean(x),min(x)*.8,sd(x)) # starting values
#out <- mh(x,v,min=min(x),cs=c(.6,.01,.01),N=50000)
##out <- out[-c(1:9000),]
#par(mfrow=c(3,1));plot(out[,1],type='l');plot(out[,2],type='l');plot(out[,3],type='l')
#apply(out,2,mean)
#v
#c(a,m,s)
#
#plot(density(x))


my.solve <- function(H,g){
  D = -H[1,1]*H[1,2]*H[1,2] + H[1,1]*H[1,1]*H[2,2] +
       H[1,3]*H[1,3]*H[2,2] + 2 * H[1,2]*H[1,3]*H[2,3] -
       H[1,1]*H[2,3]*H[2,3];

  Hig <- rep(0,3)
  Hi <- matrix(0,3,3)

  Hi[1,1] = H[1,1] * H[2,2] - H[2,3] * H[2,3];
  Hi[1,2] = H[1,3] * H[2,3] - H[1,1] * H[1,2];
  Hi[1,3] = H[1,2] * H[2,3] - H[1,3] * H[2,2];
  Hi[2,2] = H[1,1] * H[1,1] - H[1,3] * H[1,3];
  Hi[2,3] = H[1,2] * H[1,3] - H[1,1] * H[2,3];
  Hi[3,3] = H[1,1] * H[2,2] - H[1,2] * H[1,2];
  Hi[2,1] = H[1,2];
  Hi[3,1] = H[1,3];
  Hi[3,2] = H[2,3];

  Hig[1] = (Hi[1,1]*g[1] + Hi[1,2]*g[2] + Hi[1,3]*g[3])/D;
  Hig[2] = (Hi[2,1]*g[1] + Hi[2,2]*g[2] + Hi[2,3]*g[3])/D;
  Hig[3] = (Hi[3,1]*g[1] + Hi[3,2]*g[2] + Hi[3,3]*g[3])/D; 
  
  Hig
}  
