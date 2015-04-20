X <- c(294,247,267,358,423,311,450,534,438,
       697,688,630,709,627,615,999,1022,1015,
       700,850,980,1025,1021,1200,1250,1500,1650)

Y<- c(30,32,37,44,47,49,56,62,68,
      78,80,84,88,97,100,109,114,117,
      106,128,130,160,97,180,112,210,135)

plot(X,Y)

lY <- log(Y)
X2 <- X^2
plot(X,lY)

mod1 <- lm(lY~X)
mod2 <- lm(lY~X+X2)

f <- function(x){
  mod2$coef[1] + mod2$coef[2]*x + mod2$coef[3]*x^2
}

abline(mod1,col='red',lwd=3)
curve(f(x), from=200,to=1600,col='blue',lwd=3,add=T)


summary(mod1)
summary(mod2)
