my.color <- function(dat,from,to) {
  if (is(dat)[1] == "function") {
    color.fn(dat,from,to)
  } else if (is(dat)[1] == "density") {
    color.den(dat,from,to)
  } else if (is(dat)[1] == "matrix") {
    color.emp(dat,from,to)
  }
}


color.den <- function(den,from,to,col="red") {
  # Colors area under a density within an interval
  # den has to be a density object
  polygon(c(from, den$x[den$x>= from & den$x <= to], to),
          c(0, den$y[den$x>=from & den$x <= to], 0),col=col)
}

color.fn <- function(f,from,to,col="red") {
  x <- seq(from,to,by=(to-from)/10000)
  polygon(c(from, x,    to),
          c(0, f(x), 0),col=col)
}


color.emp <- function(M,from,to,col="red") {
  x <- M[,1]
  y <- M[,2]
  polygon(c(from, x[x> from & x < to],    to),
          c(0, y[x> from & x < to], 0),col=col)
}
