#1
a.image <- function(Q,color=paste0("gray",0:100),...) {
  image(t(apply(Q,2,rev)),yaxt="n",xaxt="n",col=color,...)
}

#2: Change this matrix. Put an awesome 0-1 pattern.
M1 <- matrix(c(0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,1,1,1,
               0,0,0,0,1,0,
               0,0,0,0,1,0),6,6,byrow=T)

# Use this to view your pattern
a.image(M1)

#3:

genData <- function(M1) {
  X <- matrix(0,10,36)

  m <- c(M1)
  X <- t(apply(X,1,function(x) m+rnorm(36,0,.5)))

  X
}

X <- genData(M1)

name <- readline("Enter your name here: ")
write.table(X,sep="\t",file=paste0(name,".dat"),quote=F,col.names=F,row.names=F)

# Email the .dat to luiarthur@gmail.com
# Thanks! You're freaking awesome! %) - a Picasso smiley.
