a.image <- function(Q,color=rev(heat.colors(100)),#paste0("gray",100:0),
                    numbers=F,num.cex=1,label.axis=F,
                    numcolor="black",axis.num=T,...) {

  image(t(apply(Q,2,rev)),yaxt="n",xaxt="n",col=color,...)
  
  ec <- apply(Q,2,sum) 
  er <- apply(Q,1,sum)
  seq1 <- seq(0,1,len=length(ec))
  seq2 <- seq(1,0,len=length(er))
  seq4 <- seq(0,1,len=nrow(Q))
  if (label.axis) axis(4,at=seq4,lab=1:n,las=2,...)

  if (axis.num) {
    axis(1,at=seq1,lab=ec)
    axis(2,at=seq2,lab=er,las=2,...)
  }

  if (numbers) {
    xx <- rep(1:ncol(Q),each=nrow(Q))
    yy <- rep(1:nrow(Q),ncol(Q))
    text(seq1[xx],seq2[yy],c(Q),col=numcolor,font=2,cex=num.cex)
    #print(t(Q)[xx,yy])
    #for (x in 1:ncol(Q)) {
    #  for (y in 1:nrow(Q)) {
    #    text(seq1[x],seq2[y],t(Q)[x,y],col=numcolor,font=2,cex=num.cex)
    #  }
    #}  
  }
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

plot.post <- function(x,main=NULL,hpd=T,color="cornflowerblue",cex.l=1,trace=T,stay=F,tck.dig=4,its=length(x),...) {
  mn.x <- round(mean(x),5)
  v.x <- round(sd(x),3)
  den <- density(x)
  rng <- c(min(den$y),max(den$y))

  diff <- rng[2]-rng[1]
  main <- ifelse(is.null(main),"Posterior Distribution",
                         paste("Posterior Distribution \n for",main))
  if (hpd) {
  } else {
  }

  rng.x <- range(den$x)
  x.diff <- rng.x[2] - rng.x[1]

  if (hpd) {
    hpd <- get.hpd(x)

    plot(density(x),col=color,ylim=c(rng[1],rng[2]+diff*.3),lwd=3,
         main=main,xaxt="n")

    color.den(den,rng.x[1],rng.x[2],col.den=color,col.area=color,add=T)
    color.den(den,hpd[1],hpd[2],col.den=col.mult(color),
              col.area=col.mult(color),add=T) 
    lines(c(mn.x,mn.x),c(0,bound(mn.x,den,ret=F)),lwd=2,col="red")

    axis(1,at=c(hpd,mn.x),labels=round(c(hpd,mn.x),tck.dig),las=0,...)
    legend("topleft",legend=c(paste("Mean =",mn.x),
                              paste("Std. Dev. =",v.x),
                              paste("Low HPD =",round(hpd[1],4)),
                              paste("Upp HPD =",round(hpd[2],4)),
                              paste("Iterations =",its)),
                              bty="n",cex=cex.l)
  } else {
    plot(density(x),col=color,ylim=c(rng[1],rng[2]+diff*.3),lwd=3,main=main)
    color.den(den,rng.x[1],rng.x[2],col.den=color,col.area=color,add=T)
    lines(c(mn.x,mn.x),c(0,bound(mn.x,den,ret=F)),lwd=2,col="red")
    legend("topleft",legend=c(paste("Mean =",mn.x),
                              paste("Std. Dev. =",v.x),
                              paste("Iterations =",length(x))),
                              bty="n",cex=cex.l)
  }

  mfg <- par()$mfg

  if (trace) {
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
  }

  if (!(stay)) {
    row.num <- mfg[1]
    col.num <- mfg[2]
    last.row <- mfg[3]
    last.col <- mfg[4]

    if (col.num < last.col) {
      mfg[2] <- mfg[2] + 1
    } else {
      if (row.num < last.row) {
        mfg[1] <- mfg[1] + 1
      } else {
        mfg[1] <- 1
      }
      mfg[2] <- 1
    }
  }

  par(mfg=mfg)
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
  mar <- x1 <- x2 <- y1 <- y2 <- NULL
  s <- par("usr")
  if (is.numeric(coords)) {
    x1 <- coords[1]; x2 <- coords[2]
    y1 <- coords[3]; y2 <- coords[4]
  } else if (coords=="topright"){
    x1 <- s[1] + (s[2]-s[1]) * (1-scale)
    x2 <- s[2] - (s[2]-s[1]) * .01
    y1 <- s[3] + (s[4]-s[3]) * (1-scale)
    y2 <- s[4]
    mar <- c(.1,.1,1,.1)
  } else if (coords=="bottomright") {
    x1 <- s[1] + (s[2]-s[1]) * (1-scale)
    x2 <- s[2] - (s[2]-s[1]) * .01
    y1 <- s[3] + (s[4]-s[3]) *.05
    y2 <- y1 + (s[4]-s[3]) * (scale)
    mar <- c(1,.1,1,.1)
  } else if (coords=="topleft") {
    x1 <- s[1] + (s[2]-s[1]) * .05
    x2 <- x1 + (s[2]-s[1]) * (scale)
    y1 <- s[3] + (s[4]-s[3]) * (1-scale)
    y2 <- s[4]
    mar <- c(.1,1,1,.1)
  }else if (coords=="bottomleft") {
    x1 <- s[1] + (s[2]-s[1]) * .05
    x2 <- x1 + (s[2]-s[1]) * (scale)
    y1 <- s[3] + (s[4]-s[3]) *.05
    y2 <- y1 + (s[4]-s[3]) * (scale)
    mar <- c(1,1,1,.1)
  }
  opts <- par(no.readonly=T)
    par(fig = c(grconvertX(c(x1,x2),from="user",to="ndc"),
                grconvertY(c(y1,y2),from="user",to="ndc")),
        mar = mar, new = TRUE)
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

est.Z <- function(Zs,p=.5) {
  B <- length(Zs)
  EZ <- sum.matrices(Zs) / B
  EZ <- ifelse(EZ>p,1,0)
  col0.ind <- which(apply(EZ,2,function(x) sum(x)==0))
  EZ <- as.matrix(EZ[,-col0.ind])
  EZ
}

Rapply <- function(L,f) { # L is a list, f is a function to apply to L[[x]]
                          # L apply takes a list, applies a function, 
                          # and rbinds it. Assumes output is vector.
  n <- length(L)
  out <- apply(matrix(1:n),1,function(i) f(L[[i]]))
  t(out)
}

clust.Z <- function(z) {
  z <- as.matrix(z)
  v.z <- apply(z,1,function(x) toString(x))
  uniq.vz <- unique(v.z)
  clust.num <- apply(matrix(v.z),1,function(x) which(uniq.vz %in% x))
  k <- length(uniq.vz)
  n <- nrow(z)
  z.out <- matrix(0,n,k)
  for (i in 1:n) {
    z.out[i,clust.num[i]] <- 1
  }
  z.out
}


det <- function(x,log=F) {
  out <- 0
  if (!log) {
    out <- det(x)
  } else {
    out <- unlist(determinant(x,log=T))[1]
  }
  out
}

plot.contour <- function(M,...) {
  library(MASS) # filled.contour, kde2d
  J <- kde2d(M[,1],M[,2])
  contour(J,...)
}

plot.posts <- function(M,names=rep(NULL,ncol(as.matrix(M))),digits=4,cex.legend=.7,
                       keep.par=F,tck.dig=4,cex.a=1/ncol(M),its=nrow(M),...) {
  M <- as.matrix(M)
  k <- ncol(M)
  corrs <- cor(M)
  set <- par(no.readonly=T)
  par(mfrow=c(k,k))
    for (i in 1:k) {
      if (i>1) {
        for (j in 1:(i-1)) { 
          plot(1, type="n", axes=F, xlab="", ylab="",
               main=paste0("Corr (",names[i],",",names[j],")")) # empty plot
          r <- round(corrs[i,j],digits)
          cex.cor <- max(.8/strwidth(format(r)) * abs(r),1)
          text(1,labels=r,cex=cex.cor)
          #legend("center",legend=corrs[i,j],
          #       title=paste0("Corr (",names[i],",",names[j],")"))
        }  
      }
      
      plot.post(M[,i],cex.l=cex.legend,main=names[i],tck.dig=tck.dig,cex.axis=cex.a,its=its,...)

      if (i<k) {
        for (j in (i+1):k) {
          plot(M[,c(j,i)],type="l",col="gray85",xlab=names[j],ylab=names[i],
               main=paste("Trace & Contour \n",names[i],"vs",names[j]))
          plot.contour(M[,c(j,i)],add=T)
        }
      }  
    }
  if (!(keep.par)) par(set)
}


