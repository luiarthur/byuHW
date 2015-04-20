core <- c(1294,1279,1274,1264,1263,1254,1251,
          1251,1248,1240,1232,1220,1218,1210)
peri <- c(1284,1272,1256,1254,1242,1274,1264,1256,1250)
n <- length(core)
m <- length(peri)
w <- c(core,peri)
Sp <- (sum((core-mean(core))^2) + sum((peri-mean(peri))^2) )/(n+m-2)

t.stat <- ( mean(core) - mean(peri) ) / sqrt(Sp * (1/n + 1/m) )

2*pt(t.stat,n+m-2)
