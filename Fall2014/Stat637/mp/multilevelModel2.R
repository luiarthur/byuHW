# Paper:
#http://www.stat.columbia.edu/~gelman/research/published/multi2.pdf
# Data taken from:
#http://www.unc.edu/courses/2007spring/enst/562/001/docs/assignments/assign10.htm
#http://www.radon.com/radon/radon_levels.html

# (radon) activity in this dataset was measured in pCi/L. 
# What is an acceptable level of radon gas? (Taken from Radon.com)
#   Radon Act 51 passed by Congress set the natural outdoor level of radon gas
#   (0.4 pCi/L) as the target radon level for indoor radon levels. Unfortunately
#   two-thirds of all homes exceed this level. The US EPA was tasked with setting
#   practical guidelines and recommendations for the nation. To this end, the US
#   EPA has set an action level of 4 pCi/L. At or above this level of radon, the
#   EPA recommends you take corrective measures to reduce your exposure to radon
#   gas. This does not imply that a level below 4.0 pCi/L is considered
#   acceptable, as stated in the BEIR VI study. It is estimated that a reduction
#   of radon levels to below 2 pCi/L nationwide would likely reduce the yearly
#   lung cancer deaths attributed to radon by 50%. However, even with an action
#   level of 2.0 pCi/L, the cancer risk presented by radon gas is still hundreds
#   of times greater than the risks allowed for carcinogens in our food and
#   water.

source("countdown.R")
source("plotpost.R")

radon <- read.csv("radon.txt",header=T)
city <- read.csv("cty.txt",header=T)

rad <- radon[which(radon$state=="MN"),c("activity","basement","county")]
rad <- rad[which(sapply(rad$basement,as.character) > ""),]
rad[which(rad$activity == 0),"activity"] <- .05
ctyInfo <- city[which(city$st=="MN"),c("cty","Uppm")]
dat <- merge(rad,ctyInfo,by.x="county",by.y="cty")

N <- nrow(dat)
J <- length(unique(dat$county))
y <- log(dat$activity)
x <- ifelse(dat$basement=="Y",1,0)
counties <- unique(dat$county)

G <- matrix(0,N,J)
u <- log(unique(dat$Uppm))
Y <- as.list(1:J)
X <- as.list(1:J)
for (j in 1:J) {
  Y[[j]] <- y[which(dat$county==counties[j])]
  X[[j]] <- x[which(dat$county==counties[j])]
}


mh <- function(B=1e3,csa=rep(.1,J),csb=.1,csg0=.1,csg1=.1,cssy=.1,cssa=.01,
               sigb2=100,ubg0=100,ubg1=100,ubsy2=100,ubsa2=100) {
  # Likelihoods:
  # yij ~ N(aj+b*xij, sy2)  ......(1)
    ll1 <- function(j,aj,b,sy2) sum(dnorm(Y[[j]],aj+b*X[[j]],sqrt(sy2),log=T))
  #  aj ~ N(g0+g1*uj, sa2) ......(2)
    ll2 <- function(aj,g0,g1,uj,sa2) dnorm(aj,g0+g1*uj,sqrt(sa2),log=T)
  # y ~ prod_{j=1:J}{i=1:N}{N(aj+b*xij, sy2)}  ......(3)
    ll3 <- function(a,b,sy2) {
      out <- 0
      for (j in 1:J) out <- out + sum(dnorm(Y[[j]],a[j]+b*X[[j]],sqrt(sy2),log=T))
      out
    }
  #  a ~ prod_{j=1:J}{N(g0+g1*uj, sa2)} ......(4)
    ll4 <- function(a,g0,g1,sa2) sum(dnorm(a,g0+g1*u,sqrt(sa2),log=T))

  #Priors:
  # aj ~ N(g0+g1*uj,sa2)
    laj <- ll2
  # b ~ N(0,100)
    lb <- function(b) dnorm(b,0,sqrt(sigb2),log=T)
  # g0 ~ N(0,100)
    lg0 <- function(g0) dnorm(g0,0,10,log=T) #0
  # g1 ~ N(0,100)
    lg1 <- function(g1) dnorm(g1,0,10,log=T) #0
  # sy2 ~ N(2,1)
    lsy2 <- function(sy2) dnorm(sy2,2,log=T)#0
  # sa2 ~ N(2,1)
    lsa2 <- function(sa2) dnorm(sa2,2,log=T)#0

  a <- matrix(0,B,J)
  b <- rep(0,B)
  g0 <- rep(1,B)
  g1 <- rep(1,B)
  sy2 <- rep(1,B)
  sa2 <- rep(1,B)

  acc.a <- rep(0,J)
  acc.b <- 0
  acc.g0 <- 0
  acc.g1 <- 0
  acc.sy2 <- 0
  acc.sa2 <- 0
  
  for (z in 2:B) {
    ot <- Sys.time()

    a[z,] <- a[z-1,]
    b[z] <- b[z-1]
    g0[z] <- g0[z-1]
    g1[z] <- g1[z-1]
    sy2[z] <- sy2[z-1]
    sa2[z] <- sa2[z-1]

    # Update a:
    for (j in 1:J) {
      cand <- rnorm(1,a[z,j],sqrt(csa[j]))
      lr <- ll1(j,cand,  b[z],sy2[z]) + laj(cand,  g0[z],g1[z],u[j],sa2[z]) -
           (ll1(j,a[z,j],b[z],sy2[z]) + laj(a[z,j],g0[z],g1[z],u[j],sa2[z]))
      if (lr > log(runif(1))) {
        a[z,j] <- cand
        acc.a[j] <- acc.a[j]+1
      }
    } # End of J Loops

    # Update b:
    cand <- rnorm(1,b[z],csb)
    lr <- ll3(a[z,],cand,sy2[z]) + lb(cand) - 
       (ll3(a[z,],b[z],sy2[z]) + lb(b[z]))
    if (lr > log(runif(1))) {
      b[z] <- cand
      acc.b <- acc.b+1
    }

    # Update sy2:
    cand <- rnorm(1,sy2[z],cssy)
    if (cand>0){
      lr <- ll3(a[z,],b[z],cand) + lsy2(cand) - 
           (ll3(a[z,],b[z],sy2[z]) + lsy2(sy2[z]))
      if (lr > log(runif(1))) {
        sy2[z] <- cand
        acc.sy2 <- acc.sy2+1
      }
    }

    # Update sa2:
    cand <- rnorm(1,sa2[z],cssa)
    if (cand>0){
      lr <- ll4(a[z,],g0[z],g1[z],cand) + lsa2(cand) - 
           (ll4(a[z,],g0[z],g1[z],sa2[z]) + lsa2(sa2[z]))
      if (lr > log(runif(1))) {
        sa2[z] <- cand
        acc.sa2 <- acc.sa2+1
      }
    }

    # Update g0:
    cand <- rnorm(1,g0[z],csg0)
    lr <- ll4(a[z,],cand,g1[z],sa2[z]) + lg0(cand) - 
         (ll4(a[z,],g0[z],g1[z],sa2[z]) + lg0(g0[z]))
    if (lr > log(runif(1))) {
      g0[z] <- cand
      acc.g0 <- acc.g0+1
    }

    # Update g1:
    cand <- rnorm(1,g1[z],csg1)
    lr <- ll4(a[z,],g0[z],cand,sa2[z]) + lg1(cand) - 
         (ll4(a[z,],g0[z],g1[z],sa2[z]) + lg1(g1[z]))
    if (lr > log(runif(1))) {
      g1[z] <- cand
      acc.g1 <- acc.g1+1
    }
    count.down(ot,z,B)
  }# End of B Loops
  
  out <- list("a"=a,"b"=b,"g0"=g0,"g1"=g1,"sy2"=sy2,"sa2"=sa2,
              "acc.a"=acc.a/B,"acc.b"=acc.b/B,"acc.g0"=acc.g0/B,"acc.g1"=acc.g1/B,
              "acc.sy2"=acc.sy2/B,"acc.sa2"=acc.sa2/B)
  out
}

out <- mh(B=1e5)

pdf("latex/images/apost.pdf")
  plot.posts(out$a[,c(1,50,85)],names=c("a1","a50","a85"),color="pink",cex.legend=.5)
dev.off()

pdf("latex/images/hyperPost.pdf",width=19,height=13)
  plot.posts(cbind(out$b,out$g0,out$g1,out$sy2,out$sa2),
             names=c("b","g0","g1","sy2","sa2"))
dev.off()

out$acc.a
out$acc.b
out$acc.g0
out$acc.g1
out$acc.sy2
out$acc.sa2

a.h <- apply(out$a,2,mean)
b.h <- mean(out$b)
g0.h <- mean(out$g0)
g1.h <- mean(out$g1)
y1 <- apply(as.matrix(a.h),1,function(a) a+b.h)
y0 <- apply(as.matrix(a.h),1,function(a) a)

cty <- sapply(counties,as.character)
y.m <- cbind(as.factor(cty),y1,y0,u)
y.m <- y.m[order(y.m[,2],y.m[,3]),]

pdf("latex/images/ym.pdf")
  pmar <- par("mar")
  par(mar=c(5.2,4.1,1,1))
  plot(y.m[,2],type="l",col="red",lwd=3,ylab="log radon",xaxt="n",xlab="",ylim=range(c(y1,y0,u)))
  lines(y.m[,3],col="blue",lwd=3)
  lines(y.m[,4],col="grey30",lwd=3)
  axis(1,at=1:J,lab=cty[y.m[,1]],las=2,cex.axis=.5)
  #axis(1,at=1:J,lab=paste0(cty[y.m[,1]],": ",round(y.m[,4],3)),las=2,cex.axis=.5)
  legend("topleft",legend=c("Basement","No Basement","log uranium"),col=c("red","blue","grey30"),bty="n",lwd=3)
  axis(2,at=log(4),las=2,lab="log(4)",cex.axis=.6)
  #axis(2,at=c(log(2),log(4)),las=2,lab=c("log(2)","log(4)"),cex.axis=.6)
  abline(h=log(4),col="yellow",lwd=2)
  #abline(h=log(2),col="green",lwd=2)
  par(mar=pmar)
dev.off()

pdf("latex/images/au.pdf")
  plot(u,a.h,main="aj vs. log(Uranium)",pch=20,col="grey30",
       ylab="Regression Intercept (aj)",xlab="County-level Log Uranium Meaure (uj)")
  abline(g0.h,g1.h,lwd=2,col="purple")
  legend("topleft",bty="n",
         legend=c(paste("Intercept:",round(g0.h,4)),
                  paste("Slope:",round(g1.h,4))))
dev.off()

clay <-  dat[which(dat$county=="CLAY"),]
lyon <-  dat[which(dat$county=="LYON"),]

mean(clay$act)
mean(lyon$act)

mean(clay$act>4)
mean(lyon$act>4)
