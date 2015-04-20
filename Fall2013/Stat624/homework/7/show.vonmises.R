show.vonmises <- function(n,mu,vu,kappa1,kappa2,lambda,show.contours=TRUE){
  source('../7/dvonmises.R',chdir=TRUE)
  source('../7/rejection-sampler.R',chdir=TRUE)
  source('../7/rvonmises.R',chdir=TRUE)
  
  plot(0,cex=0, xlim=c(-pi,pi), ylim=c(-pi,pi))

  if (show.contours){
    x <- seq(-pi,pi,by=.1)

    M <- matrix(0,length(x),length(x))
    for (i in 1:length(x)){
      for (j in 1:length(x)){
        X <- c(x[i],x[j])
        M[i,j] <- dvonmises(X,mu,vu,kappa1,kappa2,lambda,log=T,F)
      }
    }
    contour(x,x,M) 
  }

  samp <- rvonmises(n,mu,vu,kappa1,kappa2,lambda)
  points(samp,col='red',cex=.5)

  #samp2 <- read.table('../8/results.txt',header=F);
  #points(samp2,pch=20,col='blue')
}

#show.vonmises(1000,mu=pi,vu=pi/2,kappa1=20,kappa2=10,lambda=28,show.contours=T)
#Cdata <- read.table("../8/temp.txt")
#points(Cdata,col='seagreen3')

#system.time(show.vonmises(2000,mu=pi,vu=pi/2,kappa1=20,kappa2=10,lambda=28,show.contours=T))
