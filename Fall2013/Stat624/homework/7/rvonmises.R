rvonmises <- function(n,mu,vu,kappa1,kappa2,lambda){
  source('../7/dvonmises.R',chdir=T)
  source('../7/rejection-sampler.R',chdir=T)
  
  #phi,psi
  mode <- c(0,0)
  if (kappa1*kappa2 > lambda^2) { # unimodal
    mode <- c(mu,vu)
  } else { # bimodal
      phi0 <- acos(kappa1/abs(lambda) *
              sqrt((lambda^2+kappa2^2)/(lambda^2+kappa1^2))) + mu
      psi0 <- acos(kappa2/abs(lambda) *
              sqrt((lambda^2+kappa1^2)/(lambda^2+kappa2^2))) + vu
      if (lambda>0){
        mode <- c(phi0,psi0)
      } else {
        mode <- c(-phi0,psi0)
      }
      if (dvonmises( mode,mu,vu,kappa1,kappa2,lambda,use.norm=F) < 
          dvonmises(-mode,mu,vu,kappa1,kappa2,lambda,use.norm=F)){
        mode <- -mode
      }  
  }  

  log.f <- function(x){
      dvonmises(x,mu,vu,kappa1,kappa2,lambda,log=TRUE,use.norm=FALSE)
  }

  log.g <- function(x){
        log(1/(4*pi^2))
  }

  alpha <- exp(log.g()-log.f(mode))

  sampler <- function(){
    runif(2,-pi,pi)
  }

  samp <- matrix(0,n,2)
  samp <- rejection.sampler(log.f,log.g,sampler,alpha,n)
  
  return(samp)
}

# TEST CODE:
#n <- 1000; mu <- pi/2; vu <- pi; kappa1 <-10; kappa2 <- 10; lambda <- 24
#samp <- rvonmises(n,mu,vu,kappa1,kappa2,lambda)
#plot(samp)
