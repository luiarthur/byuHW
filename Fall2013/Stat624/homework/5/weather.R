# %s/dim(wData)\[1\]/n/gc
# READ and ORGANIZE DATA
  library(sqldf)
  wData <- read.csv('ny12.dat', header=T)
  n <- dim(wData)[1]
  # 274 was obtained by inspection of wData
  wData$Month <- substr(as.Date(wData$MST),6,7) 
  attach(wData)
  cn <- colnames(wData)

  moBegin <- NULL; moBegin[1] <- 1
  moEnd   <- NULL; moEnd[12]  <- n
  j <- 2; k <- 1
  for (i in 2: (length(Month)-1)){
    if(Month[i] != Month[i-1]) {
      moBegin[j] <- i 
      j <- j + 1
    }
  }
  moEnd[1:11] <- moBegin[2:12] - 1
  mo <- cbind(moBegin,moEnd)
  rm(i,j,k,moBegin,moEnd)


# Monthly Precip:
  sumPrec <- 
  sqldf('
          select Month, sum(PrecipitationIn) as sumPrec from wData
            group by Month
        ')
  sumPrec <- sumPrec[,2] * 100

  nData <- read.csv('normal.dat',header=T)
  nData$nMonthDay <- substr(as.Date(nData$nMST),6,10)
  subData <- 
    sqldf('
           select nMonthDay, avg(nMaxTemperatureF) as avMax, 
                  avg(nMinTemperatureF) as avMin from nData
           group by nMonthDay 
          ')

  nData$nMonth <- substr(as.Date(nData$nMST),6,7)
  nData$nYear <- substr(as.Date(nData$nMST),1,4)

  precData <- 
    sqldf('
            select nMonth, nYear, sum(nPrecipitationIn) as sumPrec from nData
              group by nYear, nMonth
          ')
  mPrecData <- 
    sqldf('
            select avg(sumPrec) as avSumPrec from precData group by nMonth
          ')
  mt <- mean(Mean.TemperatureF)
  nmt <- mean(nData$nMeanTemperatureF,na.rm=T)
  
# Daily:
plots <- function(){
    pdf('./plot.pdf')
    par(mfrow=c(3,1),xaxs='i',yaxs='i', mar=c(1,5,1,4), font=2)
    layout(matrix(c(1,1,2,3),4,1))

    plot(Max.TemperatureF,type='l',
         xlim=c(0,n),ylim=c(0,110),axes=F,
         ylab='Temperature', xlab='Month',
         main = 'Provo City\'s Weather for 2012')
    lines(Min.TemperatureF)
    polygon(c(1:n,n:1),
            c(Max.TemperatureF,rev(Min.TemperatureF)),
            col='pink')
    
    abline(h = (0:11*10), v = c(mo[,1],n))
    lines(subData$avMax,col='red')
    lines(subData$avMin,col='blue')
    centerMonth <- apply(mo[1:12,],1,mean)
    xy <- matrix(c(centerMonth,rep(105,12)),12,2)
    text(xy,c('Jan','Feb','Mar','Apr',
               'May','Jun','Jul','Aug',
               'Sep','Oct','Nov','Dec'))
    legend(5,97,legend=c('Annual Temperature','2012: 52.8','Normal: 50.5'))
    # axis(1), axis(2), axis(3), axis(4): 
    # Puts axis on bottom, left, top, right
    axis(2, at=((0:10)*10), labels=((0:10)*10), las=2)
    # NEED DAILY NORMAL MAX/MIN
    x <- which.max(Max.TemperatureF)
    y <- Max.TemperatureF[x]
    val <- bquote(.(y))
    text(x,y,paste('HIGH:',val),col='red',font=2)
    x <- which.min(Min.TemperatureF)
    y <- Min.TemperatureF[x]
    val <- bquote(.(y))
    text(x,y,paste('LOW:',val),col='blue',font=2)
    
    text(300,80,'Normal HIGH',col='red')
    text(230,40,'Normal LOW',col='blue')

    barplot(t(data.matrix(cbind(sumPrec[,2],mPrecData))), beside=T, xlab='Month', ylab='Precipitation', ylim=c(0,2))
    abline(h=0:4*20)
    legend(15.3,1.8,legend=c('Precipitation (In)',
                           'Total in 2012: 3.29   ',
                           'Normal: 73'))
  # Yes, Jan one year was an outlier? 

  # Humidity:
    plot(Mean.Humidity,type='l', axes=F, ylab='Noon Relative Humidity (%)',
         ylim=c(0,100))
    abline(h=(0:4)*25, v=c(mo[,1],n))
    axis(2, at=(0:4*25), labels=(0:4*25), las=2)
    dev.off()
}
