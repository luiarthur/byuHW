system('./clean')
system('./compile')
source('../7/dvonmises.R',chdir=TRUE)
source('../7/rvonmises.R',chdir=TRUE)
source('../7/show.vonmises.R',chdir=TRUE)
source('rvonmises2.R')
dyn.load("rvonmises2.so")

getInfo <- function(draws,Cdraws,mu,nu,k1,k2,l,show.plot=TRUE){

  if (show.plot){
    show.vonmises(0,mu,nu,k1,k2,l)
    title(main='Contour Plot with Random Draws \n Red -- R \n Cyan -- C') 
  }

  getRInfo <- function(n){
    R.time <- system.time( Rdata <- rvonmises(n,mu,nu,k1,k2,l) )
    if (show.plot) {points(Rdata,col='red',cex=.5)}
    R.time
      #list(time=R.time,draws=Rdata)
  }

  getCInfo <- function(n){
    C.time <- system.time(  Cdata <- rvonmises2(n,mu,nu,k1,k2,l) )
    if (show.plot) {points(Cdata,col='cyan',cex=.5)}
    C.time
    #list(time=C.time,draws=Cdata)
  }
  
  R.time <- rbind(t(draws),apply(draws,1,getRInfo))
  C.time <- rbind(t(draws),apply(draws,1,getCInfo))
  C.2    <- rbind(t(Cdraws),apply(Cdraws,1,getCInfo))

  rownames(R.time)[1] <- 'Number.of.Draws'
  rownames(C.time)[1] <- 'Number.of.Draws'
  rownames(C.2)[1] <-    'Number.of.Draws'

  list(R.time=R.time, C.time=C.time, C.2=C.2)
  
}

draws <- as.matrix(seq(100,3000,length=30))
Cdraws <- as.matrix(seq(1000,1000000,length=31))

#draws <- as.matrix(1:30 * 30); Cdraws <- as.matrix(1:30 * 100)
#draws <- as.matrix(100:101); Cdraws <- as.matrix(1000:1001)

# ELAPSED TIME MATRIX
# This function, should take about 3 to 5 minutes to run.
M <- getInfo(draws,Cdraws,pi,pi/2,20,10,28,show.plot=F)

r.div.c <- M$R.time[4,] / M$C.time[4,]
ratio  <- lm(r.div.c ~ M$R.time[1,])

pdf("analysis.pdf")
  par(mfrow=c(4,1))
  plot(M$R.time[1,],M$R.time[4,],main='R',xlab='Number of Draws',ylab='Time(s)',type='o')
  plot(M$C.time[1,],M$C.time[4,],main='C',xlab='Number of Draws',ylab='Time(s)',type='o')
  plot(M$R.time[1,],M$R.time[4,] / M$C.time[4,],
       main='R/C',xlab='Number of Draws',ylab='Time(s)',type='o')
  abline(ratio)
  plot(M$C.2[1,],M$C.2[4,],main='C.2',xlab='Larger Number of Draws',ylab='Time(s)',type='o')
dev.off()

#system('./clean')

#rmod <- lm(M$R.time[4,] ~ M$R.time[1,])
#cmod <- lm(M$C.time[4,] ~ M$C.time[1,])

#ratio <- rmod$coefficients[2] / cmod$coefficients[2]


# To Jared, the man.
# Here is my analysis on the speed of my rvonmises (R) and rvonmises2 (C) functions.
# From the plot generated, "analysis.pdf", we can see that the time taken to 
# execute both rvonmises (R) and rvonmises (C) increases linearly as the number
# of draws increases. The speed ratio R code and C code converges to is 130. That is,
# in the time that the R code produces ONE draw, the C code has produced 130 draws!
# BOOYA!
# Also, I know that the plot is reproduceable, but I thought you might not  
# run my code as it takes a few minutes... so I gave you the plot anyways...
