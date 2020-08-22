x<-1:5
y<-6:10
z<-11:15
ls()
sum(x,y,z)
plot(x,y)
rm(f, m, z)
BOD
date()
log(x)
log2(x)
exp(x)
ln(x)
x1=c(1,3,5,7,9)
gender=c("male", "female")
2:7
seq(from=1, to=50, by=10)
seq(from=0, to=50, by=10)
rep(1,times=10)
x2=rep(1,times=10)
x3=rep("Naimul", times=5)
rep(1:3, times=4)
rep(seq(from=2,to=7,by=0.25), times=3)
#if two vectors are of the same length, we may add/subtract/mult/div
#corresponding elements
x1[3]
x1[-3]
x[1:3]
x[c(1,5)]
x[-c(1,3)]
x[x>3]
matrix(c(1:9), nrow = 3, byrow = TRUE)
matrix(c(1:9), nrow = 3, byrow = FALSE)
mat<-matrix(c(1:9), nrow = 3, byrow = TRUE)
mat
mat[1,2]
mat[c(1,3),2]
mat[2,]
mat*100
help("read.csv")
data1 <- read.csv(file.choose(), header = T)
data2 <- read.table(file.choose(), header = T, sep=",")
#data1 & data2 are same(basically for comma seperated csv files) 
data3 <- read.delim(file.choose(), header = T)
data4 <- read.table(file.choose(), header = T, sep="\t")
#data3 & data4 are same(basically for tab delimineted txt files)
data3 <- read.table(file="E:/hs0.csv", header = T, sep=",")
rm(data3)
dim(data1)
head(data1)
#first 6 data list
tail(data1)
#last 6 data list
names(data1)
data1$read
mean(data1$read)
attach(data1)
#Should not have any collapsing element on the environment
read
mean(read)
detach(data1) 
#by attach(starts manipulation) and detach(finishes manipulation) command we can easily manipulate the variables as separate arrays
class(gender)
levels(prgtype)
summary(data1)
z<- c(0,1,1,1,0,0,0,0,0,0)
z<- as.factor(z)
#to make that variable a factor variable
science[20:38]
science[37]
data1[20:38, ]
mean(write[gender==1])
fdata <-data1[gender==0, ]
mdata <-data1[gender==1, ]
gender<-as.factor(gender)
summary(gender)
MaleOver15W <-data1[gender==1 & write>55, ]
#Multiple condition adding...
MaleOver15W[1:4, ]
set.seed(1234)
n<- 10000
x<-runif(n)
y<-runif(n)
mpi10<-4*sum(x^2+y^2<=1)/10000
ave(mpi:mpi9) -> averaged_pi
you<-10000
for(i in 1:you){
n<- 10000
x<- runif(n)
y<- runif(n)
print(i<-c(nai[i]<-4*sum(x^2+y^2<=1)/10000))}





