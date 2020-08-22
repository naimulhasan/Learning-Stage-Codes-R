x<-c(42,65,75,78,87,42,45,68,72,90,19,24,80,81,81,36,54,69,77,84,42,51,57,59,78,51,74,75,78,132,60,60,72,95,138,18,20,27,42,60,15,30,39,62,84,69,109,113,118,153,64,90,93,109,112,61,78,94,109,136)
length(x)
data<-matrix(x,nrow=12, byrow = T)
data
rownames(data)=1:12
data<-as.data.frame(data)
colnames(data)<-c("s1","s2","s3","s4","s5")
data$SMean<-apply(data,1,mean)
apply(data,2,mean)
data$Range<-apply(data[,-6],1,function(x)diff(range(x)))
data
xbar<-mean(data$SMean)
Rbar<-mean(data$Range)
A2=0.577
D3=0
D4=2.115
ucl.xbar=xbar+A2*Rbar
lcl.xbar=xbar-A2*Rbar
cl.xbar=xbar
plot(1:12,data$SMean,xlab = "Sample Number", ylab = "Sample Mean", main = "Control chart for Mean")
abline(h=cl.xbar)
abline(h=ucl.xbar,lty=3)
abline(h=lcl.xbar,lty=3)
abline(h=ucl.xbar,lty="dotted")
abline(h=lcl.xbar,lty="dotted")
legend(4,45,"Lower Control Line",text.width = 4,text.font = 8)
legend(4,80,"Central Line")
legend(4,115,"Upper Control Line")
qnorm(.98)
x<-c(.98,.92,.84,.75,.50,.20,.05)
v_Z<-function(x)list(qnorm(x))
v_Z(x)
y<-c(v_Z(x))
xy<-as.numeric(as.character(unlist(y)))
list(50+10*(xy))
year<-c(1964:1984)
sales<-c(80,84,80,88,98,92,84,88,80,100,84,96,92,104,116,112,102,114,108,126,121)
length(sales)

e29<-data.frame(year, sales)
e29_1<-e29[1:7,]
e29_2<-e29[8:14,]
e29_3<-e29[15:21,]
s1<-sum(e29_1$sales)
s2<-sum(e29_2$sales)
s3<-sum(e29_3$sales)
a=84.81
b=84.57
c=1.22
t<-c(1:21)
Ut<-(a+b*c^t)
library(ggplot2)
e<-data.frame(year,sales,Ut)

ggplot(e,aes(x=c(1:nrow(e)), y = Ut)) + geom_line()

dt<-c(12.16,12.63,13.46,14.12,14.94,15.34,15.65,17.04,17.62,18.04,18.44,18.85,18.77,19.11,19.91,20.38,20.44,20.20,20.44)
sum(log10(dt))
m<-c(.8372,.8324,.8318,.8344,.8346,.8332,.8340,.8344,.8308)
pnorm(-4.7982)

(.84-.8345)/0.002889
pnorm(1.903773)-pnorm(-5.019038)

p<-c(.01,.02,.03,.04,.05,.06,.08,.10,.30,.8)
np<-225*p

list(1-ppois(q = 13, lambda = np))
pp<-ppois(13,np)
plot(p,pp,type = "l")

