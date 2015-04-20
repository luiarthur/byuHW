
X <- 10; M <- 0; V <- 1

Args <- c(X,M,V)
source("./random-normals") # out is the output for random.normals

#command <- paste("Rscript ./random-normals",X,M,V)
#out <- system(command,intern=T)

pdf('./box-muller.pdf')
  curve(dnorm(x,0,1),from=-5,to=5,col='red',lwd=3,
        ylim=c(0,max(dnorm(M,M,V),density(out)$y)),
        main='Random Normal Density',ylab='Density')
  lines(density(out),col='blue',lwd=3)
  legend(-5,.4, legend=c('Theoretical pdf','Simulated pdf'),
         col=c('red','blue'),lwd=3)
dev.off()
