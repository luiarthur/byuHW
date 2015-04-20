data(iris)

quant <- iris[,1:4]
overallMean <- apply(quant, 2, mean)
setosaMean <- apply(quant[which(iris$Species == 'setosa'),], 2, mean)
versicolorMean <- apply(quant[which(iris$Species == 'versicolor'),], 2, mean)
virginicaMean <- apply(quant[which(iris$Species == 'virginica'),], 2, mean)


pdf('lengthWidth.pdf')
  plot((iris$Sepal.Length), (iris$Sepal.Width),
       col = rep(c('red','blue','green'),c(50,50,50)),
       xlab = 'Sepal Length', ylab = 'Sepal Width',
       main = 'Sepal Width against Sepal Length')
  legend(6.7,4.3,legend=c('setosa','versicolor','virginica'),
         col=c('red','blue','black'), lwd=3)
dev.off()
