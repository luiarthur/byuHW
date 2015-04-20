x <- system("ls | grep dat",intern=T)
Y <- matrix(0,0,36)

label <- NULL
for (i in 1:length(x)) {
  tab <- read.table(x[i])
  label[((i-1)*10+1):(i*10)] <- substr(x[i],1,nchar(x[i])-4)
  Y <- rbind(Y,tab)
}

Y <- as.matrix(Y)
write.table(Y,sep="\t",file="Y.dat",quote=F,col.names=F,row.names=F)
system("mv Y.dat ../")
