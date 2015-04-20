source("final.R")

M <- matrix(0,4,3)

M[1,] <- mean.via.direct.sampling    (0,1,100000)
M[2,] <- mean.via.rejection.sampling (0,1,100000)
M[3,] <- mean.via.importance.sampling(0,1,500000)
M[4,] <- mean.via.mcmc               (0,1,100000)

rownames(M) <- c("Direct","Rejection","Importance","MCMC")

M
