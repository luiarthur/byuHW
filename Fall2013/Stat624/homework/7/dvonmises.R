#Evaluates the density of the sine model of the bivariate von Mises distribution
#for the supplied parameters at the angle pair given by x.
dvonmises <- function(x,mu,vu,kappa1,kappa2,lambda,log=FALSE,
                      use.normalizing.constant=TRUE){

    phi <- x[1]; psi <- x[2]
    
    C <- 1
    if (use.normalizing.constant){
      m <- 1:100
      C <- 1 / (sum(choose(2*m,m) * (lambda^2/(4*kappa1*kappa2))^m *
                  besselI(kappa1,m) * besselI(kappa2,m)) * 4 * pi^2)
    }

    out <- C * exp( kappa1*cos(phi-mu) + kappa2*cos(psi-vu) +
                    lambda*sin(phi-mu)*sin(psi-vu) )

    if (!log) {
      return(out)
    } else {
      return(log(out))
    }
}
#Test values
#x <- c(1,1); mu <- 2; vu <- 2; kappa1 <- 2; kappa2 <- 2; lambda <- 2
#dvonmises(x,mu,vu,kappa1,kappa2,lambda)

