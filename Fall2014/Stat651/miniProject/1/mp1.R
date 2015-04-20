# 1
dat1 <- read.table("binomial2.dat",header=F)
colnames(dat1) <- c("date","num.bats","num.homeruns")

n <- dat1[,2]
y <- dat1[,3]

a <- 9
b <- 1

# Prior:
curve(dbeta(x,a,b),0,1,col="red",lwd=3,ylim=c(0,23),ylab="Density",xlab=expression(theta))

# Posterior:
curve(dbeta(x,a+sum(y),b+sum(n-y)),col="blue",lwd=3,add=T)

legend("topright",legend=c("Prior","Posterior","MLE"),col=c("red","blue","orange"),lwd=3)

abline(v=sum(y)/sum(n),col="orange",lwd=3)


curve(dbeta(x,a+sum(y),b+sum(n-y)),col="blue",lwd=3,add=F,ylim=c(0,45)) # Original
curve(dbeta(x,90+sum(y),10+sum(n-y)),col="darkgreen",lwd=3,add=T)
curve(dbeta(x,100+sum(y),900+sum(n-y)),col="gold",lwd=3,add=T)
curve(dbeta(x,5+sum(y),5+sum(n-y)),col="pink",lwd=3,add=T)


