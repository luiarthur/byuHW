# Arthur Lui
# Code and Answers for HW1

library(doMC)
registerDoMC(strtoi(system("nproc",intern=T))/2)

# 1a) 
# i) Marginal Density of y:
f <- function(y){
  1/(6*sqrt(2*pi)) * {exp(-((y-66)/3)^2/2) + exp(-((y-70)/3)^2/2)}
} 

#ii) Plot for density of y :
pdf("1a.pdf")
  curve(f,53,82,main="Marginal Density of y",xname="y",col="blue",lwd=3)
dev.off()

#1b)
po <- function(y,s){
  a <- (y-66)/s
  b <- (y-70)/s
  exp(-a^2/2) / (exp(-a^2/2) + exp(-b^2/2)) 
}

#1c)
set <- 1:8
pdf("1c.pdf")
  for (i in set) {
    fun <- function(y) po(y,i)
    add <- ifelse(i==1,F,T)
    curve(fun,40,100,add=add,xlab="y",ylab="f(theta|y)",main="Theta = 1",col=i,lwd=2)
  }
  legend("topright",legend=c(paste("s=",set)),col=set,lwd=3,bty="n")
dev.off()


# 3)
one.sim <- function() {
  # Get Times when Patients enter clinic:
  in.time <- rexp(1000,1/10)
  in.time <- cumsum(in.time)
  in.time <- in.time[in.time < 420]

  l <- length(in.time)
  visit.time <- runif(l,5,20)
  wait.time <- rep(0,l)
  out.time <- rep(0,l)
  pat.doc <- rep(0,l)
  doc.free.time <- c(0,0,0)

  # Do simulation (see patients) while before 4pm:
  for (i in 1:l){
    m <- which.min(doc.free.time)
    wait.time[i] <- max(doc.free.time[m] - in.time[i], 0)
    out.time[i] <- in.time[i] + wait.time[i] + visit.time[i]
    doc.free.time[m] <- out.time[i]
  }
  
  out <- c(length(in.time),sum(wait.time>0),mean(wait.time[wait.time>0]),max(out.time))
  out[4] <- ifelse(out[4]<420, 420, out[4])
  out <- matrix(out,1,4)
  colnames(out) <- c("Num.Patients","Num.Waits","Avg.Wait (mins)","Close.Time")

  out
  #list("in.times"=in.time,"wait.times"=wait.time,"out.times"=out.time)
}

# 3a) 
A3 <- one.sim()
A3[,4] <- convert.to.time(A3[,4])
print(A3,quote=F)

# 3b)
convert.to.time <- function(x,start=9) {
  hour <- x %/% 60 + 9
  min <- x %% 60
  min <- floor(min)
  min <- ifelse(min<10,paste0(0,min),min)
  paste0(hour,":",min)
}

sims <- foreach(i=1:10000,.combine=rbind) %dopar% one.sim()
sims.int <- apply(sims,2,function(x) quantile(x,c(.1,.5,.9),na.rm=T))
sims.int[,4] <- convert.to.time(sims.int[,4])
print(sims.int,quote=F)

# Answers for Q3: ######################################################

# 3a)
# Num.Patients  Num.Waits  Avg.Wait (mins)   Close.Time
# 34            4          3.17258447900527  16:00 


# 3b)
# Percentile  Num.Patients  Num.Waits  Avg.Wait (mins)   Close.Time
# 10%         34            1          1.79487111102708  16:00     
# 50%         42            5          3.91915571136854  16:06     
# 90%         51            12         6.68772664459924  16:15 

