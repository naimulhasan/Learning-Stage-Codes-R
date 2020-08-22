nm<-read.csv("G:/Project4/POSTNEO.csv",sep = ",")
nm
nmM<-read.csv("G:/Project4/POSTNEONATALMALE.csv",sep = ",")
nmF<-read.csv("G:/Project4/POSTNEONATALFEMALE.csv",sep = ",")
par(mfrow=c(3,1))
boxplot(nmF$NATIONAL,horizontal = T)
boxplot(nmF$RURAL,horizontal = T)
boxplot(nmF$URBAN,horizontal = T)

nm1<-nm[1:7,]
nm2<-nm[8:14,]
par(mfrow=c(1,1))
boxplot(nm1$NATIONAL,horizontal = T)
boxplot(nm2$NATIONAL,horizontal = T)
mean(nm2$NATIONAL)
acf(nm2$NATIONAL)
head(nm2$NATIONAL)
head(nm1$NATIONAL)
var.test(nm1$NATIONAL,nm2$NATIONAL)
t.test(nm1$NATIONAL,nm2$NATIONAL,var.equal = F)
var.test(nm1$RURAL,nm2$RURAL)
t.test(nm1$RURAL,nm2$RURAL,var.equal = T)
var.test(nm1$URBAN,nm2$URBAN)
t.test(nm1$URBAN,nm2$URBAN,var.equal = T)

plot(nm1$RURAL,nm1$URBAN)
cor.test(nm2$RURAL,nm2$URBAN)



t.test(nmM$NATIONAL,nmF$NATIONAL)
t.test(nmM$RURAL,nmF$RURAL)
t.test(nmM$URBAN,nmF$URBAN)














library(readxl)
install.packages("xlsx")
library(xlsx)
TF<-read_xlsx("G:/Project2/TFNew.xlsx",col_names = T)
TFT<-read_xlsx("G:/Project2/TFTotal.xlsx",col_names = T)
rtf<-lm(TFT$ForeignEx~TFT$Total)
summary(rtf)
TF1<-TFT[1:15,]
TF2<-TFT[16:29,]
rt1<-lm(TF1$ForeignEx~TF1$Total)
rt2<-lm(TF2$ForeignEx~TF2$Total)
summary(rt1)
summary(rt2)
coefficients(rt1)
coefficients(rt2)
TFT$PerVisEarned<-TFT$ForeignEx/TFT$Total

plot(TF1$ForeignEx, TF1$Total,main= "Corelation Between ",xlab="Foreign Ex",ylab="Total Tourist", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
plot(TF2$ForeignEx, TF2$Total,main= "Corelation Between ",xlab="Foreign Ex",ylab="Total Tourist", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
plot(TFT$ForeignEx, TFT$Total,main= "Corelation Between ",xlab="Foreign Ex",ylab="Total Tourist", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
cor(TFT$ForeignEx, TFT$Total)
cor(TF1$ForeignEx, TF1$Total)
cor(TF2$ForeignEx, TF2$Total)
cor.test(TFT$ForeignEx, TFT$Total)
cor.test(TF1$ForeignEx, TF1$Total)
cor.test(TF2$ForeignEx, TF2$Total)
View(TF1)
CV<-function(mean,sd)
{
  (sd/mean)*100
}
heights <- c(151, 160, 162, 155, 154, 168, 153, 158, 157, 150, 167)
weights <- c(61, 69, 73, 65, 64, 78, 63, 68, 67, 60, 77)

library(raster)
cv(heights)
cv(weights)
{
  
  plot(TF1$ForeignEx, TF1$Total, xlab = "Foreign Exchange)",
       ylab = "Total Tourists", col = "blue", lwd = 2)
  
  abline(rt1, col = "red")
}
summary(rtf)
rtf
{
  
  plot(TFT$ForeignEx, TFT$Total, xlab = "Foreign Exchange)",
       ylab = "Total Tourists", col = "blue", lwd = 2)
  
  abline(rtf, col = "red")
}

coefT<-coefficients(rtf)
estForeign<-coefT[1]+coefT[2]*TFT$Total
estForeign
cv(estForeign)
cv(TFT$ForeignEx)
cv(TFT$ForeignEx)*(sqrt(1-0.04856))

cv(TF1$ForeignEx)
cv(TF2$ForeignEx)
cv(TF1$ForeignEx)*(sqrt(1-0.901))
cv(TF2$ForeignEx)*(sqrt(1-0.5181))

t.test(TF0$Total,TF2$Total,var.equal=T,alternative = "less")
t.test(TF0$Total,TF2$Total,var.equal=T)
t.test(NT0$PerVisEarned,NT2$PerVisEarned,var.equal=T)
t.test(NT0$PerVisEarned,NT2$PerVisEarned,var.equal=T, alternative = "less")
TF0<-TF1[2:15,]
t.test(TF0$ForeignEx,TF2$ForeignEx,var.equal=T,alternative = "less")
t.test(TF0$ForeignEx,TF2$ForeignEx,var.equal=T)

tft<-lm(TFT$Total~TFT$Year)
summary(tft)
#R-sq: 0.1892
cv(TFT$Total)*(sqrt(1-0.1892))
tfft<-lm(TFT$ForeignEx~TFT$Year)
summary(tfft)
#R-sq: 0.8653
cv(TFT$ForeignEx)*(sqrt(1-0.8653))

tf1<-lm(TF1$Total~TF1$Year)
summary(tf1)
#R-sq: 0.7865
cv(TF1$Total)*(sqrt(1-0.7865))

tf2<-lm(TF2$Total~TF2$Year)
summary(tf2)
#R-sq: 0.5429
cv(TF2$Total)*(sqrt(1-0.5429))

ft1<-lm(TF1$ForeignEx~TF1$Year)
summary(ft1)
#R-sq: 0.7817
cv(TF1$ForeignEx)*(sqrt(1-0.7817))

nn<-lm(NT$PerVisEarned~NT$Year)
summary(nn)
cv(NT$PerVisEarned)*(sqrt(1-0.6707))

nn1<-lm(NT1$PerVisEarned~NT1$Year)
summary(nn1)
cv(NT1$PerVisEarned)*(sqrt(1-0.7911))

nn2<-lm(NT2$PerVisEarned~NT2$Year)
summary(nn2)
cv(NT2$PerVisEarned)*(sqrt(1-0.7516))



ft2<-lm(TF2$ForeignEx~TF2$Year)
summary(ft2)
#R-sq: 0.8302
cv(TF2$ForeignEx)*(sqrt(1-0.8302))
mean(TF1$Total)
mean(TF2$Total)
mean(TF1$ForeignEx)
mean(TF2$ForeignEx)
mean(NT1$PerVisEarned)
mean(NT2$PerVisEarned)

library(car)
install.packages("lmtest")
library(lmtest)
nnn<-durbin.watson(nn)
nnn
nnn1<-durbin.watson(nn1)
nnn1
nnn2<-durbin.watson(nn2)
nnn2

td1<-durbin.watson(tf1)
td1
td2<-durbin.watson(tf2)
td2
fd1<-durbin.watson(ft1)
fd1
fd2<-durbin.watson(ft2)
fd2
dtf<-durbin.watson(tft)
dtf
dtff<-durbin.watson(tfft)
dtff

library(plyr)
ddply(TFT,transform,
      Growth=c(NA,exp(diff(log(TFT$Total)))-1))
install.packages("growthrate")
library(growthrate)
growth(TFT$Total)
?`growthrate-package`

attach(TFT)
NT<-data.frame(Year,Total,ForeignEx,((ForeignEx/Total)*100000))
View(NT)
detach(TFT)
fix(NT)
attach(TF0)
NT0<-data.frame(Year,Total,ForeignEx,((ForeignEx/Total)*100000))
detach(TF0)
fix(NT0)
attach(TF1)
NT1<-data.frame(Year,Total,ForeignEx,((ForeignEx/Total)*100000))
detach(TF1)
fix(NT1)
attach(TF2)
NT2<-data.frame(Year,Total,ForeignEx,((ForeignEx/Total)*100000))
detach(TF2)
fix(NT2)


growth()
ng<-lm(log(NT$PerVisEarned)~NT$Year)
summary(ng)
ng1<-lm(log(NT1$PerVisEarned)~NT1$Year)
summary(ng1)
ng2<-lm(log(NT2$PerVisEarned)~NT2$Year)
summary(ng2)

gt<-lm(log(TFT$Total)~TFT$Year)
summary(gt)
gt1<-lm(log(TF1$Total)~TF1$Year)
summary(gt1)
gt2<-lm(log(TF2$Total)~TF2$Year)
summary(gt2)
gf<-lm(log(TFT$ForeignEx)~TFT$Year)
summary(gf)
gf1<-lm(log(TF1$ForeignEx)~TF1$Year)
summary(gf1)
gf2<-lm(log(TF2$ForeignEx)~TF2$Year)
summary(gf2)

ncv<-function(x)
{
  s<-sd(x)
  m<-mean(x)
  return((s/m)*100)
}

cvt<-function(y,x)
{
  md<-lm(y~x)
  d<-durbinWatsonTest(md)$dw
  p<-summary(md)$coefficient[2,4]
  r<-summary(md)$r.squared
  n<-ncv(y)
  t<-n*sqrt(1-summary(md)$r.squared)
  return(list("P-Value"=p,"R square"= r, "Coefficient of Variation"=n, "Coefficient of Variation Around Trend Line"=t,"Durbin Watson Statistic"=d))
}
cvt(TFT$ForeignEx,TFT$Year)

library(car)

x<-192757.1-141754.4
(x/141754.4)+1
y<-10-5
(y/5)+1
