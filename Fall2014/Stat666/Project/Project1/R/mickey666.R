dat = read.table("oliver3b.txt", head=TRUE)
dat = as.matrix(dat)[-c(3,40,176),]
n = nrow(dat)
p = ncol(dat)

mw.pairs <- function(x){
  par.params <- par()
  require(MASS)
  par(mfrow=rep(ncol(x), 2), mar=rep(0, 4))
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      if (i == j){
        hist(x[,i], axes = FALSE, main = "", col = "gray")
        legend("topright", legend = i, box.lty = 0, cex = 1.5)
      } else if(i > j){
        z = kde2d(x[,i], x[,j])
        plot(NA, xlim = range(z$x), ylim = range(z$y), axes = FALSE)
        .filled.contour(x=z$x, y=z$y, z=z$z, levels=seq(min(z$z), max(z$z), length=20), col=gray(seq(0.0, 1.0, length=20)))
      } else {
        plot(x[,j], x[,i], pch=20, axes = FALSE)
      }
      box()
    }
  }
  par(par.params)
}

mw.pairs(dat)

x.bar = apply(dat, 2, mean)
S = var(dat)
S.inv = solve(S)

# qq plotting
D = double(n)
for (i in 1:length(D))
    D[i] = as.vector(t(dat[i,]-x.bar) %*% S.inv %*% (dat[i,]-x.bar))

dev.off()
plot(qbeta((1:n - 0.5)/n, p/2, (n - p - 1)/2), n/(n-1)^2*sort(D),
    pch = 20)
abline(0,1)

### box-cox
# do the transformation
lam.func = function(x, lam){
    out = matrix(0, nrow(x), ncol(x))
    for (i in 1:ncol(x)){
        if (lam[i] == 0){
            out[,i] = log(x[,i])
        } else {
            out[,i] = (x[,i]^lam[i] - 1)/lam[i]
            }
        }
    return (out)
    }
# multivariate function to maximize
max.func = function(x, lam){
    y = lam.func(x, lam)
    -n/2*determinant((n-1)/n*cov(y))$modulus[1] +
        sum(apply(log(x), 2, sum) * (lam-1))
    }
# R function to do the random walk maximization
max.loop = function(x, lambda, niter = 1000, temperature = 0.99, width, print = FALSE){
    x = as.matrix(x)
    n = nrow(x)
    p = ncol(x)
    if (missing(lambda)){
        lam = runif(p, -2 ,2)
    } else {
        lam = lambda
        }
    if (missing(width))
        width = rep(1, p)
    temp = temperature
    current = max.func(x, lam)
    window = round(niter/20)
    for (i in 1:niter){
        cand.lam = lam
        for (j in 1:p){
            cand.lam[j] = runif(1, lam[j]-width[j], lam[j]+width[j])
            candidate = max.func(x, cand.lam)
            if (candidate > current){
                current = candidate
                lam[j] = cand.lam[j]
            } else {
                width[j] = width[j] * temp
                }
            }
        if (floor(i/window) == i/window && print)
            cat(i,"/",niter, "\n")
        }
    return (lam)
    }

mult.lam = max.loop(dat, print = TRUE)

max.func(dat, mult.lam)

mw.pairs(dat)
mw.pairs(lam.func(dat, mult.lam))

dev.off()
plot(dat[,4], dat[,1], pch=20)
#identify(dat[,4], dat[,1])
# observation 40 looks to be an outlier

plot(dat[,1], type='l')

new.dat = lam.func(dat, mult.lam)
x.bar = apply(new.dat, 2, mean)
S = var(new.dat)
S.inv = solve(S)

# qq plotting
D = double(n)
for (i in 1:length(D))
    D[i] = as.vector(t(new.dat[i,]-x.bar) %*% S.inv %*%
        (new.dat[i,]-x.bar))

plot(qbeta((1:n - 0.5)/n, p/2, (n - p - 1)/2), n/(n-1)^2*sort(D),
    pch = 20)
identify(qbeta((1:n - 0.5)/n, p/2, (n - p - 1)/2), n/(n-1)^2*sort(D))
abline(0,1)

### kurtosis and skewness

library(doMC)
registerDoMC(strtoi(system("nproc",intern=T))/2)

get.skew.kurt <- function(Y){
  n <- nrow(Y)
  p <- ncol(Y)
  S <- cov(Y) * (n-1) / n
  ybar <- matrix(apply(Y,2,mean),ncol=1)

  iS <- solve(S) 
  g <- matrix(0,nrow=n,ncol=n)

  # The parallel version. .97 seconds.
  doi <- function(i) {
    yi <- as.matrix(Y[i,],ncol=1)
    doj <- function(yj) t(yi-ybar) %*% iS %*% (yj-ybar)
    
    apply(Y,1,doj)
  }

  g <- foreach(j=1:n,.combine=rbind) %dopar% doi(j)

  # The sequential version. 1.5 seconds.
  #for (i in 1:n){
  #  yi <- as.matrix(Y[i,],ncol=1)
  #  for(j in 1:n){
  #    yj <- as.matrix(Y[j,],ncol=1)
  #    g[i,j] <- t(yi-ybar) %*% iS %*% (yj-ybar)
  #  }
  #}
  
  b1 <- sum(g^3) / n^2
  b2 <- sum(diag(g)^2)/n

  out <- c(b1,b2)
  names(out) <- c("skewness: b1","kurtosis: b2")
  out
}

system.time(skew.kurt <- get.skew.kurt(dat))
skew.kurt

trans.s.k <- get.skew.kurt(lam.func(dat, mult.lam))

# Next step. What is the transformed data? And what is the skewness and kurtosis?
# After removing outliers, what is the skewness?

z1 <- function(b1,Y) {
  n <- nrow(Y)
  p <- ncol(Y)
  ((p+1)*(n+1)*(n+3)*b1) / (6*((n+1)*(p+1)-6))
}

z1.score <- z1(trans.s.k[1],new.dat)
pchisq(z1.score,p*(p+1)*(p+2)/6,lower.tail=F)

z2 <- function(b2,Y) {
  n <- nrow(Y)
  p <- ncol(Y)
  (b2-p*(p+2)*(n+p+1)/n) / sqrt(8*p*(p+2)/(n-1))
}  
z2.score <- z2(trans.s.k[2],new.dat)
pnorm(z2.score,lower.tail=F)

