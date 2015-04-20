#1
# The simToy function takes in a probability vector, 
# and computes the mean number of purchases required
# to "Collect'em all".

returnAnswer1 <- function(ans){
  paste('Mean Number of Purchases: ', ans)
}

simToy <- function(prob=rep(1/4,4), N=10^4) {

  numOfToys <- length(prob)
  toy <- 1:numOfToys
  V <- NULL

  for (i in 1:N) {
    count <- rep(0, numOfToys)
    while ( !all(count >= rep(1,numOfToys)) ){
      newToy <- sample(toy,1,prob=prob)
      count[newToy] <- count[newToy] + 1
    }
    V[i] <- sum(count)
  }

  V
}
#1a Each Toy is Equally Likely.
simA <- simToy(rep(.25,4))
ans1a <- mean(simA)

#1b Probability of Drawing each Toy is NOT Equally Likely.
simB <- simToy(c(.1,.25,.25,.4))
ans1b <- mean(simB)

#2
returnAnswer2 <- function(ans){
  paste('Proportion of Customres that will purchase 14 or more boxes: ', ans)
}

ans2a <- mean(simA>14)
ans2b <- mean(simB>14)

#3
CI <- function(ans){
  paste('We are 95 percent confident that the true mean is contained in the interval: ',ans[1],' and ',ans[2])
}

ciA <- mean(simA) + c(-1,1)*qnorm(1-0.05/2)*sd(simA)/sqrt(length(simA))
ciB <- mean(simB) + c(-1,1)*qnorm(1-0.05/2)*sd(simB)/sqrt(length(simB))

#ANSWERS:

returnAnswer1(ans1a); returnAnswer1(ans1b); returnAnswer2(ans2a); returnAnswer2(ans2b)
