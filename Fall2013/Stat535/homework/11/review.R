apts <- read.table('http://grimshawville.byu.edu/azapts.txt',header=T,na.strings='.')

psr <- apts2
psr <- apts[,c(1,4,10)]
psr <- psr[order(psr[,3]),]
psr <- subset(psr,!is.na(psr[,1]))
#head(psr)

X0 <- ( psr[,3] == "CentPhoenix" )
X1 <- ( psr[,3] == "CentPhoenix" ) * psr[,2]

X2 <- ( psr[,3] == "Chandler" )
X3 <- ( psr[,3] == "Chandler" ) * psr[,2]

X4 <- ( psr[,3] == "Mesa" )
X5 <- ( psr[,3] == "Mesa" ) * psr[,2]

X6 <- ( psr[,3] == "NPhoenix" )
X7 <- ( psr[,3] == "NPhoenix" ) * psr[,2]

X8 <- ( psr[,3] == "NWCities" )
X9 <- ( psr[,3] == "NWCities" ) * psr[,2]

X10<- ( psr[,3] == "Scottsdale" )
X11<- ( psr[,3] == "Scottsdale" ) * psr[,2]


X12<- ( psr[,3] == "SWPhoenix" )
X13<- ( psr[,3] == "SWPhoenix" ) * psr[,2]

X14<- ( psr[,3] == "Tempe" )
X15<- ( psr[,3] == "Tempe" ) * psr[,2]

X <- cbind(X0,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15)

Y <- psr[,1]

betaHat <- solve( t(X) %*% X ) %*% t(X) %*% Y

oneBH <- lm(Y[1:sum(X0)]~X1[1:sum(X0)])
plot(X1[1:sum(X0)],Y[1:sum(X0)],cex=.1,col='red')
abline(oneBH)

points(X3[1281:1356],Y[1281:1356],cex=.1,col='blue')
abline(betaHat[3,1],betaHat[4,1])
