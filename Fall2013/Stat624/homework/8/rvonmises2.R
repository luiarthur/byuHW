if (!is.loaded("rvonmises2.so"))
  dyn.load("rvonmises2.so")

rvonmises2 <- function(n,mu,nu,k1,k2,l){
  #if (!is.loaded("rvonmises2.so")){
    #options('scipen'=10)
    #system('./clean')
    #system('./compile')
    #dyn.load("rvonmises2.so")
  #}  

  #command <- paste('./rvonmises2',n,mu,nu,k1,k2,l,'>|', data)
  #system(command)
  #Cdata <- read.table(data)
  
  PRINT <- 0 # 0 - Don't PRINT; 1 - PRINT
  phi <- rep(0,n)
  psi <- rep(0,n)
  Cdata <- .C("rvonmises",as.integer(n),as.double(mu),as.double(nu),
                          as.double(k1),as.double(k2),as.double(l),
                          phi=as.double(phi),psi=as.double(psi),
                          as.integer(PRINT))

  cbind(Cdata$phi, Cdata$psi)                        
}

#rvonmises2(10,pi,pi/2,10,10,28)
