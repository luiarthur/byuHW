# 2i) Full, linear model
set.seed(5)
n <- 50
x <- rnorm(n,4,1)
b <- 4
y <- rnorm(n,x*b,1)

plot(x,y,cex=2)
mod.full <- glm(y~x)
mod.int <- glm(y~1)

abline(mod.full,col='red',lwd=3)
abline(mod.int,col='blue',lwd=3)

bh <- mod.full$coef
X <- cbind(1,x)
D.1 <- sum((y - X%*%bh)^2)
mod.full
D.1
pchisq(D.1,n-1)

# 2i) Intercept, linear model
bh.int <- mod.int$coef
D.0 <- sum((y-bh.int)^2)
mod.int
D.0 
pchisq(D.0,n-1)

# 2ii) Intercept, normal
norm.int <- glm(x~1)
bh.2 <- norm.int$coef
D.0.norm <- sum((x-bh.2)^2)
D.0.norm
norm.int
pchisq(D.0.norm,n-1)

# 2ii) Full, normal
xx <- seq(n,1)
norm.full <- glm(x~xx)
bh.3 <- norm.full$coef
D.1.norm <- sum((x-cbind(1,xx)%*%bh.3)^2)
norm.full
D.1.norm
pchisq(D.1.norm,n-1)

# 2b: How do the 4 models compare to the corresponding saturated models?
p1 <- pchisq(D.0,n-1,lower=F) # significantly different from saturated model
p2 <- pchisq(D.1,n-2,lower=F) # not significantly different from sat. model
p3 <- pchisq(D.0.norm,n-1,lower=F) # not significantly different from sat. model 
p4 <- pchisq(D.1.norm,n-2,lower=F) # not significantly different from sat. model

# 2c:
d.D <- D.0 - D.1
d.D.norm <- D.0.norm - D.1.norm

pchisq(d.D,1,lower=F) # Reject: slope is better
pchisq(d.D.norm,1,lower=F) # Fail to Reject: slope does not improve model

