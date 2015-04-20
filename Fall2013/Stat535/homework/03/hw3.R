#################################################################
# Bivariate Normal Distribution                                 #
biv.norm<-function(x,y,mu.x,mu.y,sigma.x,sigma.y,rho){          #
   (1/(2 * pi * sigma.x * sigma.y * sqrt(1-rho^2))) *           #
   exp( -(1/(2*(1-rho^2))) * ( ((x-mu.x)/sigma.x)^2 +           #
         ((y-mu.y)/sigma.y)^2 - (2 * rho * ((x-mu.x)/sigma.x) * #
         ((y-mu.y)/sigma.y) ) ) )                               #
}                                                               #
#################################################################

# GENERATE DATA:

  #install.packages("mvtnorm")
  library("mvtnorm")

  #n <- 10^4
  #rho <- 1/2
  #sigX <- sqrt(2); sigY <- sqrt(2)
  #muX <- 5; muY <- 5
  #cov <- rho * sigX * sigY

  
  #X <- rmvnorm( n, c(muX,muY), matrix( c(sigX^2,cov,cov,sigY^2), 2, 2) )
  X <- read.table('hw3.dat', header=F,sep="\t")
  x <- X[,1]; y <- X[,2] 
  n <- dim(X)[1]

pdf('plot.pdf')

  plot(x,y, pch=20,cex=.25)
  title(main='E[Y|X] and E[X|Y]')
  legend(-1,11,legend=c('E[Y|X]','E[X|Y]','95% Region'),
         col=c('red','blue','green'),lwd=2)


# 95% Region
  xyz <- cbind(X,1:n,rep(0,n))
  colnames(xyz) <- c('x','y','ind','probs')
  mx <- mean(x); my <- mean(y);
  sx <- sd(x); sy <- sd(y)
  rho <- cov(x,y) / (sd(x)*sd(y))

  xyz[,4] <- biv.norm(x,y,mx,my,sx,sy,rho)

  xyz <- as.data.frame(xyz)
  xyz <- xyz[order(xyz$probs),]
  xyz <- xyz[-c(1:(.05*n)),]

  points(xyz[,1],xyz[,2],col='green',pch=20, cex=.25)

# Lines
  #mod <- lm(y~x)
  #abline(mod,col='yellow')
  # E[Y|X]
  slp <- cov(x,y)/var(x)
  int <- mean(y) - slp*mean(x)
  abline(int,slp,col='red',lwd=2)
  # E[X|Y]
  slp <- cov(y,x)/var(y)
  int <- mean(x) - slp*mean(y)
  abline(-int/slp,1/slp,col='blue',lwd=2)

dev.off()

