f <- function(x,p0) {
  xbar <- mean(x)
  n <- length(x)
  -2 * ( sum(x) * log(p0*(1-xbar)/(xbar*(1-p0))) + n * log((1-p0)/(1-xbar)))
}

n <- 100
p0 <- .8
p <- .8

one.it <- function() {
  x <- sample(0:1,n,replace=T,prob=c(1-p,p))
  f(x,p0)
}

library(foreach)
library(doMC)
registerDoMC(16)

out <- foreach(i=1:10000,.combine=rbind) %dopar% one.it()

plot.34 <- function() {  
  hist(out,add=F,prob=T,main="Simulation",col="red")
  curve(dchisq(x,1),from=0,to=10,add=T,col="blue",lwd=3)
  legend("topright",legend=c(paste("n = ",n),paste("p0 =",p0),paste("p = ",p)))
}

pdf("plot.pdf"); plot.34(); dev.off()
