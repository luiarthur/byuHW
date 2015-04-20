a.image <- function(Q,color=paste0("gray",0:100),...) {
  image(t(apply(Q,2,rev)),yaxt="n",xaxt="n",col=color,...)
}

count.down <- function(old.time,i,B) {
  prog <- round(100*i/B,4)
  new.time <- Sys.time()
  time.diff <- as.numeric(new.time-old.time)
  time.remain <- time.diff * (B-i)
  if (time.remain < 60) {
    secs <- round(time.remain)
    time.remain <- paste0(secs,"s     ")
  } else if (time.remain<3600) {
    mins <- round(time.remain%/%60)
    secs <- round(time.remain%%60)
    time.remain <- paste0(mins,"m ",secs,"s        ")
  } else {
    hrs <- round(time.remain%/%3600)
    mins <- round((time.remain%%3600) %/% 60)
    time.remain <- paste0(hrs,"h ",mins,"m         ")
  }
  cat(paste0("\rProgress: ",prog,"%. Time Remaining: ",time.remain," "))
  if (i==B) cat("100%\n")
}

# Example
#B <- 5000
#for(i in 1:B) {
#  old.time <- Sys.time()
#  Sys.sleep(1) # your function here
#  count.down(old.time)
#}


sum.matrices <- function(Ms,return.matrices=F) { 
# Ms is a list of matrices of different lengths
# return.matrices is a boolean. If FALSE, function returns the sum of the matrices.
# If TRUE, function returns a list of the matrices also.

  l <- length(Ms)
  max.c <- max(unlist(lapply(Ms,ncol)))
  max.r <- max(unlist(lapply(Ms,nrow)))
  
  for (i in 1:l) {
    M <- Ms[[i]]
    
    ncol0 <- max.c - ncol(M)
    nrow0 <- max.r - nrow(M)

    if (ncol0>0) {
      col0 <- matrix(0,nrow(M),ncol0)
      M <- Ms[[i]] <- cbind(Ms[[i]],col0)
    }

    if (nrow0>0) {
      row0 <- matrix(0,nrow0,ncol(M))
      M <- Ms[[i]] <- rbind(Ms[[i]],row0)
    }
  }

  out <- Reduce("+",Ms)

  if (return.matrices) out <- list("sum"=out,"matrices"=Ms)

  out
}


# EXAMPLE: #########################################################3
#A <- matrix(1:10,nrow=2)
#B <- matrix(1:6,nrow=3)
#C <- matrix(1:6,nrow=1)
#D <- matrix(1:4)
#
#Ms <- list(A,B,C,D)
#
#
#sum.matrices(Ms)
#sum.matrices(Ms,T)
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

bound <- function(x, dens, return.x=TRUE){
  # Mickey Warner: 
  # https://github.com/mickwar/r-sandbox/blob/master/mcmc/bayes_functions.R
  # returns the x-value in dens that is closest
  # to the given x
  if (return.x)
      return(dens$x[which.min(abs(dens$x-x))])

  # returns the y-value in dens at the closest x
  return(dens$y[which.min(abs(dens$x-x))])
}

col.mult = function(col1 = 0x000000, col2 = "gray50"){
  # Mickey Warner: 
  # https://github.com/mickwar/r-sandbox/blob/master/mcmc/bayes_functions.R
  # returns the x-value in dens that is closest
  # to the given x
  if (is.character(col1))
      val1 = t(col2rgb(col1) / 255)
  if (is.numeric(col1))
      val1 = t(int2rgb(col1) / 255)
  if (is.character(col2))
      val2 = t(col2rgb(col2) / 255)
  if (is.numeric(col2))
      val2 = t(int2rgb(col2) / 255)
  rgb(val1 * val2)
}

int2rgb = function(x){
# int2rgb()
# convert an integer between 0 and 16777215 = 256^3 - 1,
# or between 0 and 0xFFFFFF
# this function is depended upon by col.mult
  # Mickey Warner: 
  # https://github.com/mickwar/r-sandbox/blob/master/mcmc/bayes_functions.R
  # returns the x-value in dens that is closest
  # to the given x
  hex = as.character(as.hexmode(x))
  hex = paste0("#", paste0(rep("0", 6-nchar(hex)), collapse=""), hex)
  col2rgb(hex)
}

plot.post <- function(x,main=NULL,hpd=T,color="cornflowerblue") {
  mn.x <- round(mean(x),5)
  v.x <- round(var(x),3)
  den <- density(x)
  rng <- c(min(den$y),max(den$y))

  diff <- rng[2]-rng[1]
  main <- ifelse(is.null(main),"Posterior Distribution",
                         paste("Posterior Distribution for",main))
  plot(density(x),col=color,ylim=c(rng[1],rng[2]+diff*.3),lwd=3,
       main=main)
  legend("topleft",legend=c(paste("Mean =",mn.x),
                            paste("Variance = ",v.x)),bty="n")
  rng.x <- range(den$x)
  x.diff <- rng.x[2] - rng.x[1]

  opts <- par(no.readonly=T)
    left <- rng.x[1] + x.diff*2/3
    right <- rng.x[2]
    par(fig = c(grconvertX(c(left,right),from="user",to="ndc"),
                grconvertY(c(rng[2],rng[2]+diff*.3),from="user",to="ndc")),
        mar = c(.1,.1,1,.1), new = TRUE)
    plot(x,type="l",col="gray30",cex.main=.5,axes=F,main="Trace Plot")
    axis(1,cex.axis=.5)
    axis(2,cex.axis=.5)
  par(opts)

  color.den(den,rng.x[1],rng.x[2],col.den=color,col.area=color,add=T)
  if (hpd) {
    hpd <- get.hpd(x)
    color.den(den,hpd[1],hpd[2],col.den=col.mult(color),
              col.area=col.mult(color),add=T)
  }

  lines(c(mn.x,mn.x),c(0,bound(mn.x,den,ret=F)),lwd=2,col="red")
  #abline(v=mn.x,col="red",lwd=2)
}

get.hpd <- function(x,a=.05,len=1e3) {
  V <- matrix(seq(0,a,length=len))
  quants <- t(apply(V,1,function(v) quantile(x,c(v,v+1-a))))
  diff <- quants[,2]-quants[,1]
  min.d <- V[which.min(diff)]
  hpd <- quantile(x,c(min.d,min.d+1-a))
  hpd
}


plot.in.plot <- function(minor.plot,coords="topright",scale=1/3) {
  # coords = x1,y1,x2,y2
  # minor.plot is a function with no parameters that plots the smaller plot
  x1 <- x2 <- y1 <- y2 <- NULL
  if (is.numeric(coords)) {
    x1 <- coords[1]; x2 <- coords[2]
    y1 <- coords[3]; y2 <- coords[4]
  } else if (coords=="topright"){
    s <- par("usr")
    x1 <- s[1] + (s[2]-s[1]) * (1-scale)
    x2 <- s[2] - (s[2]-s[1]) * .01
    y1 <- s[3] + (s[4]-s[3]) * (1-scale)
    y2 <- s[4]
  } else if (coords=="bottomright") {
    s <- par("usr")
    x1 <- s[1] + (s[2]-s[1]) * (1-scale)
    x2 <- s[2] - (s[2]-s[1]) * .01
    y1 <- s[3] + (s[4]-s[3]) * .1
    y2 <- s[3] + (s[4]-s[3]) * (scale)
  }
  opts <- par(no.readonly=T)
    par(fig = c(grconvertX(c(x1,x2),from="user",to="ndc"),
                grconvertY(c(y1,y2),from="user",to="ndc")),
        mar = c(.1,.1,1,.1), new = TRUE)
    minor.plot()
    #axis(1,cex.axis=.5)
    #axis(2,cex.axis=.5)
  par(opts)
}

#x <- rnorm(10000)
#plot(density(x),ylim=c(0,.5))
#
#minor <- function() {
#  plot(x,type="l",axes=F,main="Trace",cex.main=.8) 
#  axis(1,cex.axis=.5)
#  axis(2,cex.axis=.5)
#}
#
##plotinplot(minor,c(1,4,.4,.5))
#plot.in.plot(minor,"topright")
#plot.in.plot(minor,"bottomright")

