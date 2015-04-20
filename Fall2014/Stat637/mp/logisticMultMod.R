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
thresh <- log(4) 
y <- log(dat$activity)
z <- ifelse(y>thresh,1,0)
x <- ifelse(dat$basement=="Y",1,0)
counties <- unique(dat$county)
u <- log(unique(dat$Uppm))

Y <- as.list(1:J)
Z <- as.list(1:J)
X <- as.list(1:J)
for (j in 1:J) {
  Y[[j]] <- y[which(dat$county==counties[j])]
  Z[[j]] <- z[which(dat$county==counties[j])]
  X[[j]] <- x[which(dat$county==counties[j])]
}


# Logistic Veresion
mh <- function(B=1e3,csa=rep(2,J),csb=.3,csg0=.1,csg1=.3,cssa=.05,sigb2=100) {
  # Likelihoods:
  # zij ~ Bern(pij)  ......(1)
    ll1 <- function(j,zj,aj,b) {
      pj <- exp(aj+b*X[[j]]) / (1+exp(aj+b*X[[j]]))
      sum(dbinom(zj,1,pj,log=T))
    }  
  #  aj ~ N(g0+g1*uj, sa2) ......(2)
    ll2 <- function(aj,g0,g1,uj,sa2) dnorm(aj,g0+g1*uj,sqrt(sa2),log=T)
  # y ~ prod_{j=1:J}{i=1:N}{N(aj+b*xij, sy2)}  ......(3)
    ll3 <- function(a,b) {
      out <- 0
      for (j in 1:J) {
        pj <- exp(a[j]+b*X[[j]]) / (1+exp(a[j]+b*X[[j]]))
        out <- out + sum(dbinom(Z[[j]],1,pj,log=T))
      }  
      out
    }
  #  a ~ prod_{j=1:J}{N(g0+g1*uj, sa2)} ......(4)
    ll4 <- function(a,g0,g1,sa2) sum(dnorm(a,g0+g1*u,sqrt(sa2),log=T))

  #Priors:
  # aj ~ N(g0+g1*uj,sa2)
    laj <- ll2
  # b ~ N(0,100)
    lb <- function(b) dnorm(b,0,sqrt(sigb2),log=T)
  # g0 ~ U(-100,100)
    lg0 <- function(g0) dnorm(g0,0,10,log=T) #0
  # g1 ~ U(-100,100)
    lg1 <- function(g1) dnorm(g1,0,10,log=T) #0
  # sy2 ~ U(0,100)
    lsy2 <- function(sy2) dnorm(sy2,2,log=T)#0
  # sa2 ~ U(0,100)
    lsa2 <- function(sa2) dnorm(sa2,2,log=T)#0

  a <- matrix(0,B,J)
  b <- rep(0,B)
  g0 <- rep(1,B)
  g1 <- rep(1,B)
  sa2 <- rep(1,B)

  acc.a <- rep(0,J)
  acc.b <- 0
  acc.g0 <- 0
  acc.g1 <- 0
  acc.sa2 <- 0
  
  for (l in 2:B) {
    ot <- Sys.time()

    a[l,] <- a[l-1,]
    b[l] <- b[l-1]
    g0[l] <- g0[l-1]
    g1[l] <- g1[l-1]
    sa2[l] <- sa2[l-1]

    # Update a:
    for (j in 1:J) {
      cand <- rnorm(1,a[l,j],sqrt(csa[j]))
      lr <- ll1(j,Z[[j]],cand,  b[l]) + laj(cand,  g0[l],g1[l],u[j],sa2[l]) -
           (ll1(j,Z[[j]],a[l,j],b[l]) + laj(a[l,j],g0[l],g1[l],u[j],sa2[l]))
      if (lr > log(runif(1))) {
        a[l,j] <- cand
        acc.a[j] <- acc.a[j]+1
      }
    } # End of J Loops

    # Update b:
    cand <- rnorm(1,b[l],csb)
    lr <- ll3(a[l,],cand) + lb(cand) - 
       (ll3(a[l,],b[l]) + lb(b[l]))
    if (lr > log(runif(1))) {
      b[l] <- cand
      acc.b <- acc.b+1
    }

    # Update sa2:
    cand <- rnorm(1,sa2[l],cssa)
    if (cand>0){
      lr <- ll4(a[l,],g0[l],g1[l],cand) + lsa2(cand) - 
           (ll4(a[l,],g0[l],g1[l],sa2[l]) + lsa2(sa2[l]))
      if (lr > log(runif(1))) {
        sa2[l] <- cand
        acc.sa2 <- acc.sa2+1
      }
    }

    # Update g0:
    cand <- rnorm(1,g0[l],csg0)
    lr <- ll4(a[l,],cand,g1[l],sa2[l]) + lg0(cand) - 
         (ll4(a[l,],g0[l],g1[l],sa2[l]) + lg0(g0[l]))
    if (lr > log(runif(1))) {
      g0[l] <- cand
      acc.g0 <- acc.g0+1
    }

    # Update g1:
    cand <- rnorm(1,g1[l],csg1)
    lr <- ll4(a[l,],g0[l],cand,sa2[l]) + lg1(cand) - 
         (ll4(a[l,],g0[l],g1[l],sa2[l]) + lg1(g1[l]))
    if (lr > log(runif(1))) {
      g1[l] <- cand
      acc.g1 <- acc.g1+1
    }
    count.down(ot,l,B)
  }# End of B Loops
  
  out <- list("a"=a,"b"=b,"g0"=g0,"g1"=g1,"sa2"=sa2,
              "acc.a"=acc.a/B,"acc.b"=acc.b/B,"acc.g0"=acc.g0/B,"acc.g1"=acc.g1/B,
              "acc.sa2"=acc.sa2/B)
  out
}

out <- mh(B=1e5)
pdf("latex/images/lapost.pdf")
  plot.posts(out$a[,c(1,50,85)],names=c("a1","a50","a85"),color="pink",cex.legend=.5)
dev.off()
pdf("latex/images/lhyperPost.pdf",width=19,height=13)
  plot.posts(cbind(out$b,out$g0,out$g1,out$sa2),
             names=c("b","g0","g1","sa2"))
dev.off()

out$acc.a
out$acc.b
out$acc.g0
out$acc.g1
out$acc.sa2

a.h <- apply(out$a,2,mean)
b.h <- mean(out$b)
g0.h <- mean(out$g0)
g1.h <- mean(out$g1)
p.h1 <- apply(as.matrix(a.h),1,function(a) exp(a+b.h)/(1+exp(a+b.h)))
p.h0 <- apply(as.matrix(a.h),1,function(a) exp(a)/(1+exp(a)))

cty <- sapply(counties,as.character)
p.m <- cbind(as.factor(cty),p.h1,p.h0,u)
p.m <- p.m[order(p.m[,2],p.m[,3]),]


pdf("latex/images/pm.pdf")
  pmar <- par("mar")
  par(mar=c(5.2,4.9,1,1),mfrow=c(2,1))
  #par(mar=c(1,4.9,1,1),mfrow=c(2,1))
  plot(p.m[,2],col="red",type="l",lwd=3,ylab="Probability of Radon Levels \n Exceeding 4",xaxt="n",xlab="",ylim=range(c(p.h1,p.h0)))
  lines(p.m[,3],col="blue",lwd=3)
  axis(1,at=1:J,lab=cty[p.m[,1]],las=2,cex.axis=.5)
  #axis(1,at=1:J,lab=paste0(cty[p.m[,1]],": ",round(y.m[,4],3)),las=2,cex.axis=.5)
  legend("topleft",legend=c("Basement","No Basement"),col=c("red","blue"),bty="n",lwd=3)
  plot(p.m[,4],type="l",col="grey30",lwd=3,ylab="log uranium levels (ppm)",xaxt="n",xlab="")
  axis(1,at=1:J,lab=cty[p.m[,1]],las=2,cex.axis=.5)
  par(mar=pmar,mfrow=c(1,1))
dev.off()

pdf("latex/images/lau.pdf")
  plot(u,a.h,main="aj vs. log(Uranium)",pch=20,col="grey30",
       ylab="Regression Intercept (aj)",xlab="County-level Log Uranium Meaure (uj)")
  abline(g0.h,g1.h,lwd=2,col="purple")
  legend("topleft",bty="n",
         legend=c(paste("Intercept:",round(g0.h,4)),
                  paste("Slope:",round(g1.h,4))))
dev.off()
