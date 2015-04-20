#1:
y <- 560:1000
sum(choose(1000,y) * .5^1000)

#2:
i <- 0:10
exp(-15) * sum(15^i/factorial(i))
