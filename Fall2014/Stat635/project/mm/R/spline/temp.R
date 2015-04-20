n <- 1
Z <- diag(n) %x% matrix(1,3)

t(apply(matrix(1:(n*3)),1,function(i) i*Z[i,]))
