set.seed(1)
n <- 100
X <- cbind(1,
           runif(n,-1,1),
           runif(n,-1,1),
           runif(n,-1,1),
           runif(n,-1,1))

b <- runif(5,-2,2)
p <- 1 / (1+exp(-X%*%b))
y <- rbinom(n,1,p)

x <- X[,-1]
mod <- glm(y~x,family=binomial)
cbind(mod$coef,b)

write.table(cbind(y,x),file="sim.dat",quote=F,col.names=F,row.names=F)
