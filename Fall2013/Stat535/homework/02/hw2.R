#1
news <-  c(6766, 18606, 21478, 53160, -24315, 9757)

#1a
var(news)
sum((news-mean(news))^2)/(length(news)-1)
# 629732593

#1b
n <- length(news)
mu <- mean(news)
t.stat <- qt(.975,n-1)
CI <- c(mu-t.stat * sd(news)/sqrt(n), mu+ t.stat * sd(news)/sqrt(n))

ans1b <- paste("We are 95% confident that the true mean difference lies between ", CI[1], "and ", CI[2])
cat(ans1b)

#1c
n <- length(news)
p.val <- 2 * (1 - pt( mu/(sd(news) / sqrt(n)), n-1 ))
ans1c <- paste("The p-value is: ", p.val)
cat(ans1c)
