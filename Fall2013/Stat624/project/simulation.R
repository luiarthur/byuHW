#rm(list=ls())
source('frechet.R')
dyn.load("c/frechet.so")

n <- c(500,500,1000,1000)
N <- 50
param1 <- as.double(c(8,4,1))  # a,m,s
param2 <- as.double(c(15,7,4)) # a,m,s
param <- list(param1,param2,
              param1,param2)

BigSim <- function(n=c(500,1000)[2], N=c(5,50)[1],init=param1){
  M <- matrix(0,N,3)

  M.mle  <- matrix(.C("simMLE", as.integer(N), as.integer(n), 
                   as.double(init),est=as.double(M))$est,N,3,byrow=T)

  cs    <- c(3,.1,.1)
  M.be  <- matrix(.C("simBE", as.integer(N), as.integer(n), 
                  as.double(init),as.double(cs),est=as.double(M))$est,N,3,byrow=T)

  M.me  <- matrix(.C("simME", as.integer(N), as.integer(n), 
                  as.double(init),est=as.double(M))$est,N,3,byrow=T)
  
  list(mle=M.mle,bayes=M.be,me=M.me)
}

library(foreach)
library(doMC)
registerDoMC(16)


# N =50. 
# [[1]]: n=500,  param1
# [[2]]: n=500,  param2
# [[3]]: n=1000, param1
# [[4]]: n=1000, param2
Sim <- foreach(k=1:4) %dopar% BigSim(n=n[k],N=N,init=param[[k]])

Bias.Sim <- list()
MSE.Sim  <- list()
for (i in 1:4){
  bias <- matrix(0,3,3)
  colnames(bias) <- c("a","m","s")
  rownames(bias) <- c("MLE","BE","ME")
  mse <- matrix(0,3,3)
  colnames(mse) <- c("a","m","s")
  rownames(mse) <- c("MLE","BE","ME")
  for (j in 1:3){
    bias[j,] <- apply(Sim[[i]][[j]],2,mean,na.rm=T) - param[[i]]
    mse[j,] <- apply(Sim[[i]][[j]],2,var,na.rm=T) + bias[j,]^2
  }
  Bias.Sim[[i]] <- bias
  MSE.Sim[[i]] <- mse
}
Bias.Sim
MSE.Sim

plot.bias.sim <-function(i=c(1,2)[1]){
  yl <- c(-6,20) 
  if (i==1) {
    yl <- c(-1,5)
  }

  par(mfrow=c(3,1))
  plot (n[c(i,i+2)],c(Bias.Sim[[i]][1,1],Bias.Sim[[i+2]][1,1]),type='o',col=2,
        ylab="Bias",xlab="Sample Size (n)",ylim=yl) # MLE
  lines(n[c(i,i+2)],c(Bias.Sim[[i]][2,1],Bias.Sim[[i+2]][2,1]),type='o',col=3) # BE
  lines(n[c(i,i+2)],c(Bias.Sim[[i]][3,1],Bias.Sim[[i+2]][3,1]),type='o',col=4) # ME
  legend("topright",col=2:4,lwd=3,legend=c("MLE","BE","ME"))

  plot (n[c(i,i+2)],c(Bias.Sim[[i]][1,2],Bias.Sim[[i+2]][1,2]),type='o',col=2,
        ylab="Bias",xlab="Sample Size (n)",ylim=yl) # MLE
  lines(n[c(i,i+2)],c(Bias.Sim[[i]][2,2],Bias.Sim[[i+2]][2,2]),type='o',col=3) # BE
  lines(n[c(i,i+2)],c(Bias.Sim[[i]][3,2],Bias.Sim[[i+2]][3,2]),type='o',col=4) # ME
  legend("topright",col=2:4,lwd=3,legend=c("MLE","BE","ME"))

  plot (n[c(i,i+2)],c(Bias.Sim[[i]][1,3],Bias.Sim[[i+2]][1,3]),type='o',col=2,
        ylab="Bias",xlab="Sample Size (n)",ylim=yl) # MLE
  lines(n[c(i,i+2)],c(Bias.Sim[[i]][2,3],Bias.Sim[[i+2]][2,3]),type='o',col=3) # BE
  lines(n[c(i,i+2)],c(Bias.Sim[[i]][3,3],Bias.Sim[[i+2]][3,3]),type='o',col=4) # ME
  legend("topright",col=2:4,lwd=3,legend=c("MLE","BE","ME"))
}
plot.bias.sim(1)
plot.bias.sim(2)

plot.mse.sim <-function(i=c(1,2)[1]){
  yl <- c(0,900) 
  if (i==1) {
    yl <- c(0,30)
  }

  par(mfrow=c(3,1))
  plot (n[c(i,i+2)],c(MSE.Sim[[i]][1,1],MSE.Sim[[i+2]][1,1]),type='o',col=2,
        ylab="MSE",xlab="Sample Size (n)",ylim=yl) # MLE
  lines(n[c(i,i+2)],c(MSE.Sim[[i]][2,1],MSE.Sim[[i+2]][2,1]),type='o',col=3) # BE
  lines(n[c(i,i+2)],c(MSE.Sim[[i]][3,1],MSE.Sim[[i+2]][3,1]),type='o',col=4) # ME
  legend("topright",col=2:4,lwd=3,legend=c("MLE","BE","ME"))

  plot (n[c(i,i+2)],c(MSE.Sim[[i]][1,2],MSE.Sim[[i+2]][1,2]),type='o',col=2,
        ylab="MSE",xlab="Sample Size (n)",ylim=yl) # MLE
  lines(n[c(i,i+2)],c(MSE.Sim[[i]][2,2],MSE.Sim[[i+2]][2,2]),type='o',col=3) # BE
  lines(n[c(i,i+2)],c(MSE.Sim[[i]][3,2],MSE.Sim[[i+2]][3,2]),type='o',col=4) # ME
  legend("topright",col=2:4,lwd=3,legend=c("MLE","BE","ME"))

  plot (n[c(i,i+2)],c(MSE.Sim[[i]][1,3],MSE.Sim[[i+2]][1,3]),type='o',col=2,
        ylab="MSE",xlab="Sample Size (n)",ylim=yl) # MLE
  lines(n[c(i,i+2)],c(MSE.Sim[[i]][2,3],MSE.Sim[[i+2]][2,3]),type='o',col=3) # BE
  lines(n[c(i,i+2)],c(MSE.Sim[[i]][3,3],MSE.Sim[[i+2]][3,3]),type='o',col=4) # ME
  legend("topright",col=2:4,lwd=3,legend=c("MLE","BE","ME"))
}
plot.mse.sim(1)
plot.mse.sim(2)

#EARTHQAUKE STUFF:
data(quakes)
x <- quakes$mag
init <- c(19,-1,5)
min <- min(x)
N <- 100000
n <- length(x)
outA <- rep(0,N)
outM <- rep(0,N)
outS <- rep(0,N)
#par(mfrow=c(1,1))
cs <- c(3,.10,.10)
Cdata <- .C("mig",as.double(x),as.integer(n),as.integer(N),
                  A=as.double(outA), M=as.double(outM), S=as.double(outS),
                  as.double(init),as.double(cs))
A <- Cdata$A; plot(A,type='l'); mA <- mean(A[(N*.9):N])
M <- Cdata$M; plot(M,type='l'); mM <- mean(M[(N*.9):N])
S <- Cdata$S; plot(S,type='l'); mS <- mean(S[(N*.9):N])

#plot(density(x),lwd=3,ylim=c(0,1.2))
plot.quake <- function(){  
  hist(x,prob=T,ylim=c(0,8),main="Histogram of Fiji Earthquake Magnitudes",
       xlab="Earthquake Magnitudes (Richter)")
  curve(dfrechet(x,mA,mM,mS),from=4,to=6,add=T,col='blue',lwd=3)

  meq <- .C("myEst",as.double(x),as.integer(n),est=init)$est
  lines(density(rfrechet(100000,meq[1],meq[2],meq[3])),col='green',lwd=3)

  legend("topright",legend=c("Data","Bayes Estimate","My Estimate"),
         col=c('black','blue','Green'),lwd=5)
}
plot.quake()

# Model Misspecification
# Generate Gamma(2,1), fit Frechet
misp <- rgamma(1000,3,1)
n <- length(misp)
N <- 10000
cs <- c(.5,.10,.10)
outA <- rep(0,N)
outM <- rep(0,N)
outS <- rep(0,N)
Cdata <- .C("mig",as.double(misp),as.integer(n),as.integer(N),
                  A=as.double(outA), M=as.double(outM), S=as.double(outS),
                  as.double(init),as.double(cs))
A <- Cdata$A; plot(A,type='l'); mA <- mean(A[(N*.9):N])
M <- Cdata$M; plot(M,type='l'); mM <- mean(M[(N*.9):N])
S <- Cdata$S; plot(S,type='l'); mS <- mean(S[(N*.9):N])

plot.mod.misspecify <- function() {
  hist(misp,prob=T,ylim=c(0,.3),main="Model Misspecification",xlab="Data")
  curve(dfrechet(x,mA,mM,mS),from=0,to=12,add=T,col="red",lwd=3)
  legend("topright",col=c("black","red"),lwd=3,legend=c("Gamma(3,1)","Frechet"))
}

plot.mod.misspecify()
