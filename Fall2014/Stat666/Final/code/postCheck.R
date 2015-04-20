source("rfunctions.R")

x <- rnorm(10000)
plot.post(x)

qq <- sort(apply(as.matrix(x),1,function(y) mean(y<x)))
plot(qq,type="l")

ks.test(pnorm(x),"punif")

plot(density(pnorm(x),from=0,to=1))
plot(density(pnorm(x)))

D <- apply(as.matrix(x),1,function(y) mean(y<x) - pnorm(y))

plot(sort(pnorm(x)),sort(qq))
