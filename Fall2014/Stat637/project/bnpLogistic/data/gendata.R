gendat <- function(n=100,x=seq(0.1,.99,len=n),b=c(1,-4,4),thresh=.3,sd=.1,printz=F) {
  X <- cbind(1,x,x^2)
  z <- X%*%b + rnorm(n,sd=sd)
  if (printz) print(z)
  y <- ifelse(z>thresh,1,0)
  list("x"=x,"y"=y,"z"=z,"b"=b)
}

#dat <- gendat()
#x <- dat$x
#y <- dat$y
#z <- dat$z
#
#plot(x,z,type="l",ylim=(range(c(0,1,y))))
#points(x,y)

#glm(y~x,family=binomial)
