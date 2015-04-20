# Parameters:
  # alpha > 0 (shape)
  # (Optionally, two more parameters)
  # s > 0 (scale, default: s = 1)
  # -inf < m < inf (location of min, default:m=0)

# Support:
  # x > m

# pdf:
  # alpha/x * ((x-m)/s) ^ (-1-alpha) * exp(-(x-m)/s) ^ (-alpha)

# cdf:
  # exp((-(x-m)/s)^-alpha)

# Mean:
  # m + s*Gamma(1-1/alpha), alpha > 1
  # inf, otherwise

# Median:
  # m + s * log(2) ^ (-1/alpha)

# Variance:
  # s^2 * (Gamma(1-2/alpha) - (Gamma(1-1/alpha))^2), alpha > 2
  # inf, otherwise

#rm(list=ls())

dfrechet <- function(x,a,m=0,s=1){
  # a = shape > 0; m = min; s = scale
  # x > m
  if (any(x<=m)) stop("x must be greater that m \n")
  a/s * ((x-m)/s) ^ (-1-a) * exp(-((x-m)/s) ^ (-a))
}
#Testing:
#d <- .0125 
#x <- seq(d,100000,d)
#x <- seq(d,5,d)
#sum(dfrechet(x,3)*d)

pfrechet <- function(x,a,m=0,s=1){
  if (any(x<=m)) stop("x must be greater that m \n")
  exp(-((x-m)/s)^(-a))
}

rfrechet <- function(n,a,m=0,s=1){
  if (a<=0) stop("alpha must be greater than 0\n")
  u <- runif(n)
  (-log(u))^(-1/a) * s + m
}

# Note that temp1 and temp2 are "the same"
#temp1 <- rfrechet(1000000,6,5,3)
#theoretical.stat.frechet(6,5,3)
#mean(temp1); var(temp1)
#
#temp2 <- rfrechet(1000000,6) * 3 + 5
#mean(temp2); var(temp2)

theoretical.stat.frechet <- function(a,m=0,s=1){
  mu <- ifelse(a>1, m + s*gamma(1-1/a), Inf)
  med <- m + s * log(2) ^ (-1/a)
  ss <- ifelse(a>2,s^2 * (gamma(1-2/a) - (gamma(1-1/a))^2),Inf)
  stats <- matrix(c(mu,med,ss),nrow=1)
  colnames(stats) <- c('Mean','Median','Variance')
  stats
}
#Testing: Note that when a=1, mean is supposed to be Inf.
#                             While the mean of the sample is not Inf, it is 
#                             large compared to the median.
#         Note that when a=2, var is supposed to be Inf.
#                             While the var is not Inf, it is large
#                             compared to the median.
#theoretical.stat.frechet(a)
#matrix(c(mean(samp),quantile(samp,.5),var(samp)),nrow=1)
#a <- (1:10)[5] 
#f <- function(x) dfrechet(x,a)
#samp <- rfrechet(1000000,a)
#curve(f,from=.01,to=90,lwd=3,col='red')
#lines(density(samp),col='blue',lwd=3)

#Helpful websites:
#http://www.emeraldinsight.com/journals.htm?articleid=871493&show=html#0830160102006.png
#http://en.wikipedia.org/wiki/Fr%C3%A9chet_distribution
