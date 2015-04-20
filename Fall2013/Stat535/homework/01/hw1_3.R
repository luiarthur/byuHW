> hr <- read.table("~/School/Stat535/homework/1/hw1.dat",header=T)
> hr <- hr[,9]
> mean(hr)
[1] 32.55

> (1/20)*t(rep(1,20))%*%hr
      [1,]
[,1] 32.55

Therefore, mean(hr) = (1/20)*t(rep(1,20))%*%hr = 32.55
