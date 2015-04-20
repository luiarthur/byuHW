new.urn <- function(n.balls=2,black.proportion=0.5) {
  return(c('black'=n.balls*black.proportion,
           'white'=n.balls*(1-black.proportion)))
}

proportion <- function(urn) {
  return(urn['black']/sum(urn))
}

draw.ball <- function(urn) {
  prob.black <- proportion(urn)
  drawn.ball <- c('white','black')[1*(runif(1)<prob.black)+1]
  urn[drawn.ball] <- urn[drawn.ball] + 1
  return(urn)
}

# Initial proportion of black balls
theta <- 0.5
n.balls.to.start <- 10
n.draws <- 100
precision <- 10000

doit <- function(precision) {
  draws <- numeric(precision)
  for ( i in 1:precision ) {
    urn <- new.urn(n.balls.to.start,theta)
    for ( j in 1:n.draws ) {
      urn <- draw.ball(urn)
    }
    draws[i] <- proportion(urn)
  }
  draws
}

draws <- doit(10000)

observed.test.statistic <- 0.75
p.value <- mean(draws>=observed.test.statistic)
ci <- p.value + c(-1,1)*qnorm(1-0.05/2)*sqrt(p.value*(1-p.value)/precision)
cat("p-value is ",p.value,", with a 95% confidence interval of (",ci[1],",",ci[2],")\n",sep="")

# More on graphics in Chap. 7
hist(draws,freq=FALSE,main=expression(paste("Sampling Distribution of Test Statistic under ",H[0])),xlab="Test Statistic")
lines(density(draws))
abline(v=observed.test.statistic,lwd=3)
text(observed.test.statistic+0.02,1.5,adj=0,paste("Test Stat. =",observed.test.statistic))

# 1. If alpha = .05, what is the rejection region?
     rr <- quantile(draws,c(.95,1))
     ans1 <- paste('The rejection region is the proportions between ',rr[1],' and ',rr[2],'.')
# 2. If the true proportion of black balls is .6,
#    what is the power of the test for the null hypothesis 
#    that the proportion of black balls is .5, vs. the 
#    one-sided alternative that the proportion is greater
#    than .5?
     theta <- .6
     draws2 <- doit(precision)
     power2 <- mean(draws2>=rr)
     ci2 <- power2 + c(-1,1) * qnorm(1-power2/2) * sqrt(power2*(1-power2)/precision)
     ans2 <- paste('The power if the true proportion is .6 = ', power2, 
                   ' (',ci2[1],', ',ci2[2],').',sep='')
# 3. What is the POWER if the true proportion is .7?
     theta <- .7
     draws3 <- doit(precision)
     power3 <- mean(draws3>=rr)
     ci3 <- power3 + c(-1,1) * qnorm(1-power2/2) * sqrt(power3*(1-power3)/precision)
     ans3 <- paste('The power if the true proportion is .7 = ', power3, 
                   ' (',ci3[1],', ',ci3[2],').',sep='')
## BE SURE TO ASSES MONTE CARLO ERROR!


cat('\n'); cat(ans1,'\n'); cat(ans2,'\n'); cat(ans3,'\n')
