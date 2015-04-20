my.color <- function(dat,from,to,col.den="black",col.area="red",...) {
  if (is(dat)[1] == "function") {
    color.fn(dat,from,to)
  } else if (is(dat)[1] == "density") {
    color.den(dat,from,to,col.den,col.area,...)
  } else if (is(dat)[1] == "matrix") {
    color.emp(dat,from,to)
  }
}


color.den <- function(den,from,to,col.den="black",col.area="red",add=F,...) {
  # Colors area under a density within an interval
  # den has to be a density object
  if (add) {
    #lines(den,col=col.den,...)
  } else {
    plot(den,col=col.den,...)
  }
  polygon(c(from, den$x[den$x>=from & den$x<=to], to),
          c(0, den$y[den$x>=from & den$x<=to], 0),
          col=col.area,border=col.den)
}

color.fn <- function(f,from,to) {
  x <- seq(from,to,by=(to-from)/10000)
  polygon(c(from,x,to),
          c(0,f(x),0),col="red")
}


color.emp <- function(M,from,to) {
  x <- M[,1]
  y <- M[,2]
  polygon(c(from,x[x>=from & x<= to],to),
          c(0,y[x>=from & x<=to],0),col="red")
}


#Examples: ######################################################

## color.den
#  x <- rnorm(10000)
#  denx <- density(x)
#  color.den(denx,1,2,col.area="blue")
#  color.den(denx,-1,0,col.area="red",add=T)
#
## color.fn
#  fn <- function(x) dnorm(x)
#  curve(dnorm(x),-5,5)
#  color.fn(fn,-2,5)
#
## color.emp
#  x <- seq(-5,5,length=1000)
#  y <- dnorm(x)
#  plot(x,y,type='l')
#  color.emp(cbind(x,y),-2,5)


