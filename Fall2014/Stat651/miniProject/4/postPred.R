source("../../Final/code/rfunctions.R")
source("mp4.R")

plot.post(post.pred)

qq <- apply(as.matrix(post.pred),1,function(x) mean(x<post.pred))

plot(density(qq,from=0,to=1),col=col.mult("cornflowerblue","grey80"),lwd=5)
ks.test(qq,"punif")

