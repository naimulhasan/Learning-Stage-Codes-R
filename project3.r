library(readxl)
install.packages("xlsx")
library(xlsx)
er<-read_xlsx("G:/Project3/dataset.xlsx",col_names = T)

View(er1)
names(er)<-c("year","TEx","CRem")
er1<-er[,c(1,2,4)]
save(er1,file="G:/Project3/er1.Rdata")
er<-er1
er$PExR<-er$CRem/er$TEx
er1<-er[1:21,]
er2<-er[22:42,]
mean(er1$CRem)
t.test(er1$CRem,er2$CRem,var.equal = F)
t.test(er1$TEx,er2$TEx,var.equal = F)
t.test(er1$PExR,er2$PExR,var.equal = F)
var.test(er1$CRem,er2$CRem)
var.test(er1$TEx,er2$TEx)
var.test(er1$PExR,er2$PExR)

cor.test(er$TEx,er$CRem)
cor.test(er1$TEx,er1$CRem)
cor.test(er2$TEx,er2$CRem)

r0<-lm(er$CRem~er$TEx)
r1<-lm(er1$CRem~er1$TEx)
r2<-lm(er2$CRem~er2$TEx)
summary(r0)
summary(r1)
summary(r2)
#growth rate
g0E<-lm(log(er$TEx)~er$year)
g1E<-lm(log(er1$TEx)~er1$year)
g2E<-lm(log(er2$TEx)~er2$year)
summary(g0E)
summary(g1E)
summary(g2E)
g0R<-lm(log(er$CRem)~er$year)
g1R<-lm(log(er1$CRem)~er1$year)
g2R<-lm(log(er2$CRem)~er2$year)
summary(g0R)
summary(g1R)
summary(g2R)
g0P<-lm(log(er$PExR)~er$year)
g1P<-lm(log(er1$PExR)~er1$year)
g2P<-lm(log(er2$PExR)~er2$year)
summary(g0P)
summary(g1P)
summary(g2P)
#last table
e0<-lm(er$TEx~er$year)
e1<-lm(er1$TEx~er1$year)
e2<-lm(er2$TEx~er2$year)
summary(e0)
summary(e1)
summary(e2)
cv(er$TEx)*(sqrt(1-0.7266))
cv(er1$TEx)*(sqrt(1-0.8622))
cv(er2$TEx)*(sqrt(1-0.5143))
de0<-durbin.watson(e0)
de1<-durbin.watson(e1)
de2<-durbin.watson(e2)
de0
de1
de2


r0<-lm(er$CRem~er$year)
r1<-lm(er1$CRem~er1$year)
r2<-lm(er2$CRem~er2$year)
summary(r0)
summary(r1)
summary(r2)
cv(er$CRem)*(sqrt(1-0.7121))
cv(er1$CRem)*(sqrt(1-0.9451))
cv(er2$CRem)*(sqrt(1-0.9274))

dr0<-durbin.watson(r0)
dr1<-durbin.watson(r1)
dr2<-durbin.watson(r2)
dr0
dr1
dr2


p0<-lm(er$PExR~er$year)
p1<-lm(er1$PExR~er1$year)
p2<-lm(er2$PExR~er2$year)
summary(p0)
summary(p1)
summary(p2)
cv(er$PExR)*(sqrt(1-0.6595))
cv(er1$PExR)*(sqrt(1-0.4381))
cv(er2$PExR)*(sqrt(1-0.6313))
dp0<-durbin.watson(p0)
dp1<-durbin.watson(p1)
dp2<-durbin.watson(p2)
dp0
dp1
dp2



#R-sq: 0.1892
cv(TFT$Total)*(sqrt(1-0.1892))
tfft<-lm(TFT$ForeignEx~TFT$Year)
summary(tfft)
#R-sq: 0.8653
cv(TFT$ForeignEx)*(sqrt(1-0.8653))



g1<-lm(log($PerVisEarned)~NT$Year)
summary(ng)
ng1<-lm(log(NT1$PerVisEarned)~NT1$Year)
summary(ng1)
ng2<-lm(log(NT2$PerVisEarned)~NT2$Year)
summary(ng2)



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

plot(TF1$ForeignEx, TF1$Total,main= "Corelation Between ",xlab="Foreign Ex",ylab="Total Tourist", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
plot(TF2$ForeignEx, TF2$Total,main= "Corelation Between ",xlab="Foreign Ex",ylab="Total Tourist", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
plot(TFT$ForeignEx, TFT$Total,main= "Corelation Between ",xlab="Foreign Ex",ylab="Total Tourist", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
cor(TFT$ForeignEx, TFT$Total)
cor(TF1$ForeignEx, TF1$Total)
cor(TF2$ForeignEx, TF2$Total)
cor.test(TFT$ForeignEx, TFT$Total)
cor.test(TF1$ForeignEx, TF1$Total)
cor.test(TF2$ForeignEx, TF2$Total)

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
#coef*100

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





