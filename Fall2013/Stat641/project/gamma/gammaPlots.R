to <- 15
x <- seq(from=0,to=to,by=.01)
pdf('./plot.pdf')
par(xaxs='i',yaxs='i')
curve(dgamma(x,.7,scale=2), from=0, to=to, lwd=3, col = 'red',
main = expression(GAMMA(alpha,beta)), ylab = 'Density', axes=F)
curve(dgamma(x,1,scale=2), from=0, to=to, add=T, lwd=3, col='blue')
curve(dgamma(x,2,scale=2), from=0, to=to, add=T, lwd=3, col='green')
legend(7,.45, legend=c(expression(alpha < 1), 
                      expression(alpha = 1),
                      expression(alpha > 1)),
       col=c('red','blue','green'), lwd=3, cex=2  )
axis(1, at=0, labels=0, tick=F)
abline(v=0,h=0)
axis(2, at=1/2, labels=expression(1/beta), las=2)
dev.off()

