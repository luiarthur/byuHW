## You should not modify the function definition.
## Also, you should not add any code before or after this function.
slope.simulation <- function(regressor,intercept,slope,error.distribution=c("normal","chisq","correlated")[1],nreps,alpha1,alpha2) {

  beta1Hat <- NULL
  covered <- NULL

  n <- length(regressor)
  
  for (i in 1:nreps){

    if (error.distribution == 'normal')
    {
      error <- rnorm(n) 
    } else if (error.distribution == 'chisq')
    {
      error <- rchisq(n,.5) - .5
    } else if (error.distribution == 'correlated')
    {
      mx <- mean(regressor); sx <- sd(regressor)
      error <- rnorm(n, .2 * (regressor-mx)/sx, 1)
    }

    response <- regressor * slope + intercept + error

    fm <- lm (response ~ regressor)
    beta1Hat[i] <- summary(fm)$coefficients[2,1] 
    se <- summary(fm)$coefficients[2,2]  
    CI <- beta1Hat[i] + c(1,-1) * qt(alpha1/2,n-2) * se
    covered[i] <- ifelse((CI[1] <= slope) & (slope <= CI[2]),T,F)
    
  }

  bias <- mean(beta1Hat) - slope
  coverage <- mean(covered == T)

  biasCI <- bias + c(1,-1) * qnorm(alpha2/2)*sqrt(var(beta1Hat)/nreps) 
  coverageCI <- coverage + c(1,-1) * qnorm(alpha2/2)*sqrt(coverage*(1-coverage)/nreps)

  c(bias, biasCI[1], biasCI[2], coverage, coverageCI[1], coverageCI[2])
  #NOW: I NEED TO TEX THIS UP!
}

#1) when the error terms are correlate
#2) NORM: Yes, becuase theoretically, 90% of the time, 
#   our ci should include the true slope.
#   CHISQ: No, because the 95% ci for the coverage doesn't contain .9.
#   Corr.: No, because the 95% ci for the coverage doesn't conatin .9.
#3) when n in large, an incorrect coverage is more noticable.  
