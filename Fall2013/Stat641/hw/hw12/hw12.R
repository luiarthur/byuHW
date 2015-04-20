source("../../../Stat624/project/frechet.R")

sim.Poisson <- function(mu=2,n=c(seq(3,18,length=4)),N=10000){
    par(mfrow=c(4,1))
    for (i in 1:length(n)){
      X <- matrix(rpois(N*n[i],mu),N,n[i])
      z <- ( apply(X,1,mean) - mu ) / (sqrt(mu/n[i]))
      hist(z,col=1+i,lwd=3,prob=T,xlim=c(-4,4))
      curve(dnorm(x),from=-4,to=4,lwd=3,add=T)
      legend("topright",legend=c("Normal(0,1)",paste("Poisson(",mu,")","n =",n[i])),
             col=c(1,i+1),lwd=3 )
    }
}

sim.frechet <- function(a=3,m=4,s=5,n=seq(3,999,length=4),N=10000){
    params <- theoretical.stat.frechet(a,m,s)
    mu = params[1]; ss = params[3]

    par(mfrow=c(4,1))
    for (i in 1:length(n)){
      X <- matrix(rfrechet(N*n[i],a,m,s),N,n[i])
      z <- ( apply(X,1,mean) - mu ) / sqrt(ss/n[i])
      hist(z,col=1+i,lwd=3,prob=T,breaks=20,xlim=c(-4,4),ylim=c(0,.4))
      curve(dnorm(x),from=-4,to=4,lwd=3,add=T)
      legend("topright",
             legend=c("Normal(0,1)",paste("frechet(",a,m,s,")","n =",n[i])),
             col=c(1,i+1),lwd=3)
    }
}



sim.Poisson()
sim.frechet()


