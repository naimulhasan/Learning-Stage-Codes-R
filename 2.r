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
#natural logarithm
log2(x)
exp(x)
sin((90)*(pi/180))
ar
cos((0)*(pi/180))
x1;x2
#multiple statment
round(68/3,1)

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
print(matrix(c(1:9), nrow = 3, byrow = TRUE))
matrix(c(1:9), nrow = 3, byrow = FALSE)
mat<-matrix(c(1:9), nrow = 3, byrow = TRUE)
matrix(c(1:9), nrow = 3)
x<-matrix(c(45,122,67,38), nrow = 2, byrow =T)
chisq.test(x)
x<-data.frame(c(45,122),c(67,38))

x[1,2]
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
data3 <- read.table(file="https://mirrors.netix.net/sourceforge/i/ir/irisdss/IRIS.csv", header = T, sep=",")
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
new_pi<-ave(nai)
nai<-c(1:10)
new_pi[1]
for(i in 1:you){
  n<- 10000
  x<- runif(n)
  y<- runif(n)
  i<-c(nai[i]<-4*sum(x^2+y^2<=1)/10000)}
rbinom()
rnorm(10,1,.5)
#for renaming the variable
student<-data.frame(regN,dept,age)
names(student)<-c("regN","dept","age")

barplot(percent, main="TITLE",xlab = "Gender",ylab = "%",las=1, names.arg = c("Female","Male"))

barplot(x, main="TITLE",ylab = "Gender",xlab = "%",las=1, names.arg = c("Female","Male"), horiz = T)

pie(x, main="TITLE")

boxplot(LungCap, main="TiTLE", ylab="Lung Capacity", ylim=c(0,16), las=1)

boxplot(LungCap~Gender, main="TiTLE", ylab="Lung Capacity", ylim=c(0,16), las=1)



mat1<-matrix(c(1:6),2,3)
#row,coloumn

xxy<-scan()
#for continue reading data
plot(age,regN)
summary(x)
#min, 1st quartile,median,mean,3rd quartile,max
x<-c(1:9)
y<-c(11:19)

plot(x,y)
hist(x)
hs=hist(x)
boxplot(x)
factorial(88)
gamma(5)
beta(2,3)
round(2.332356,3)
abs(18/-12)
# R is object-oriented programming language. 
# Output of a function is an object.
# It can be used as input of another 
# function.
class(b)
#datatype or type
class(log)
floor(x)
ceiling(x)
round(x)
sort(age)
sort(age,T)
#sort(age, decreasing=T)
y[3]=8
#element wise replacement in vector
length(y)
sum(x)
mean(x)
var(x)
sd(x)

data1 <- read.delim("G:/STA 334/fusion_33486_1_section2data/Pumpkin.dat", header = TRUE, sep="\t") 
data <- read.delim("G:/STA 334/fusion_33486_1_section2data/Temperature.dat", header = TRUE, sep="\t") 
install.packages("gmodels")
library(gmodels)
x<-c(1:10)
y<-c(11:20)
CrossTable(x,y)
sum(x^3)
mean((x-mean(x))^4)
age_reg=c(regN,age)
age1=(age>23)
#boolean type creates
sum(age>23)
#counting, not summation
which(age>23)
#index checking
v2 = (u > 2 & u <= 3)
v3 = (u < 2 | u == 3)
v4 = (u != 2)
u[u < 1 | u > 2] = 0
#Entering value condition wise
age[age<24]
#c1=scan(,"character") best option for character insertion
dept_factor<-factor(dept)
class(dept_factor)
# Factors
# For using in ANOVA etc., we can convert 
# character or numeric vectors into factors
c2=rep(0:2,each=10)
#for each individual repetition
ageW=ifelse(age<24,1,0)
#boolean conversion by prefferable digit
w=rep(1,length(age))
w = cut(u, breaks=c(10, 11.6, 12.5, 13), labels=c(0,1,2))
w = cut(u, c(10, 11.6, 12.5, 13), c(0,1,2))
#labeling in different points
w = as.numeric(w)-1
# Here, w has become a factor
# If you want numerical vector, use 

nn = 1
for (i in 1:100)
{
  nn = nn + 1
  if (factorial(nn) > 987654)
    break
}
nn
#use of for loop

# NA means missing values in R
# "NA" means a character (or name).
# NaN means "Not a number" (0/0).

Y1 = matrix(c(1:19, 23), 4)
Y1 = matrix(c(1:19, 23), 4,byrow = T)
#here 4 means row number
#byrow true means row wise distribution
matrix(1,10,5)
# 10x5 matrix of ones
matrix(0,10,5)
# 10x5 matrix of zeros
diag(10)
#identity matrix
diag(1:5)
#diagonal elements are 1 to 5
diag(c(2,5,8,9))
#matrix(,5,4)
t(Y1)
#transpose
Y1^-1
#inverse
W<- t(Y1)%*% Y1
#matrix multiplication
y1<-matrix(c(1:4),2,2)
W<-(y1)^-1%*%y1
A1<-matrix(1:20,5,4)
A2<-matrix(1,5,3)
print(A<-cbind(A2,A1))
#cbind is used for combining matrixes horizontally
AB<-rbind(A,B)
#rbind is used for combining matrixes vertically


x

bilas <- 50

for(i in 1:1){
if(bilas < 0) 
{print("Negative")}
else if (bilas == 50)
{print("its ok")}
else
{print("nothing")}
}
Xx = array(1:90, c(5, 6, 3))
#Higher Dimensional Arrays Like 2 or three dimensional
dim(Xx)
X~binomial(10,3)
dbinom(2,10,0.3)
#P(X=2)
pbinom(3,10,0.3)
#P(X<=3)
pbinom(2,10,0.3)
#P(X<3)
1-pbinom(4,10,0.3)
#P(X>4)
qbinom(0.5, 10, 0.3)  
# Binomial Example2 : X ~ binomial(10, 0.3)
# Find value k such that P(X <= k) = 0.5 

gender<-c(rep("Female",20),rep("Male",20))
gender<-rep(c("Female","Male"),20)
gender<-as.factor(gender)
gender2<-ordered(gender)
summary(gender)
data(iris)
dt1<-iris
head(dt1)
tail(dt1)
dt1[51:55,]
dt1[51:55,4]
dt1[51:55,3:4]
edit(dt1)
fix(dt1)
dt1<-dt1[,1:4]
summary(dt1)
m<-sapply(dt1,mean)
s<-sapply(dt1,sd)
#sapply is used for data.farmes
cv<-s/m*100
zsl<-(dt1[,1]-m[1])/s[1]
zsw<-(dt1[,2]-m[2])/s[2]
m1<-matrix(c(1:16),4,4)
dim(m1)
m1
#standardized
mean(zsl)
sd(zsl)
sd(zsw)
mean(zsl^3)
#skewness of 1st variable
mean(zsl^4)
#kurtosis
write.table(iris,"g:/iris.txt",sep=",")
datai<-read.table("g:/iris.txt",header = T,sep = ",")
fix(datai)
savehistory("location/jhk.txt")
loadhistory("location/jhk.txt")
a<-cov(dt1)*149#product matrix,except 149 covariance matrix
b<-solve(a)#inverse matrix
det(a) 
t(a)
ad<-diag(a)#trace
c<-cor(a)


a<-cov(dt1)
#positive, semi or full definite matrix, symetric or square matrix
a1<-149*a
b<-solve(a1)
det(a1)
t()
ad<-diag(a1)
#diagonal element er srt & inverse korte hobe
ad<sqrt(ad)
adhi<-1/ad
#R= D^-1/2 *S* D^-1/2
c<-cor(dt1)
e<-eigen(b)
prod(e$values)
det(b)
sum(e$values)
sum(diag(b))
dh<-sqrt(sv$d)
sv<-sh$u %*% diag(sh$d^2) %*% t(sh$v)
#Product of eigen value is the determinant of a matrix
#Sum of eigen value is the trace of a matrix
sh<-svd(b)#singular value dicomposition of not symetric. A=U' ^ V
sv#square matrix
sv-b%*%b
sr<-sh$u %*% diag(sh$d^.5) %*% t(sh$v)
sr#square root matrix
sr%*%sr
a
##REgression
sl1<-c(dt1$sl)
sw1<-c(dt1$sw)
pl1<-c(dt1$pl)
pw1<-c(dt1$pw)
Regression_iris<-lm(sl1~sw1+pl1+pw1)
Regression_iris1<-lm(sl1~sw1+pl1+pw1+pw1:sw1)

summary(Regression_iris)

pairs(dt1)
FF<-matrix(c((rep(1,each=150)),sw1,pl1,pw1),150,4)
FN<-matrix(c(sl1),150,1)
print(beta_1<-(solve(t(FF)%*%FF))%*%(t(FF)%*%FN))

FF<-matrix(c((sw1-mean(sw1)),(pl1-mean(pl1)),(pw1-mean(pw1))),150,3)
FN<-matrix(c(sl1-mean(sl1)),150,1)
print(beta_1<-(solve(t(FF)%*%FF))%*%(t(FF)%*%FN))



for(i in c(2,4,7,9))
  print(i^2)
#can insert individual values
for(i in c(2,4,7,9)){
     j<-i^2
     k<-i^3
     print(c(i,j,k))
}
for(i in c(2,4,7,9)){
  j<-i^2
  k<-i^3
  print(list(i,j,k))
}
data(iris)
iris[!complete.cases(iris)]
#For Checking missing ones
x<-rnorm(10)
y<-rnorm(10)
xy<-data.frame(x,y)
xy[5,1]<-NA
xy[!complete.cases(xy),]
meansd<-function(x)
{
  m<-mean(x)
  s<-sd(x)
  #print(list(m,s))
  r<-list("Mean"=m,"SD"=s)
  #return(r) can be used to display
  #or, print(r)
}

print(meansd(dt1$sl))

fix(meansd)

meansd(iris[,1])
setwd("location")
getwd()
installed.packages()
install.packages("epiR")
install.packages("psych")
install.packages("rattle")
library(psych)
library(rattle)
help(package=epiR)
install.packages()
?rexp()
rexp(25)
data1<-data.frame(iris)
data1
data1<-data1[which(data1$Species=="setosa"|data1$Species=="versicolor"),]
#specific collection
t.test(data1$Sepal.Length~data1$Species)
#for unequal variance t-test
t.test(data1$Sepal.Length~data1$Species,var.equal=T)
t.test(data1$Sepal.Length~data1$Species,var.equal=T, alternative="less")
t.test(data1$Sepal.Length~data1$Species,var.equal=T,alternative="greater")
#for equal variance t-test
t.test(data1$Sepal.Length)
#single mean test by default 0
t.test(data1$Sepal.Length,mu=5)
#single mean test by mean 5
t.test(Sepal.Length,mu=5,data=data1)
t.test(data1$Sepal.Length,data1$Sepal.Width,paired=T)
#paired t-test
t.test(data1$Sepal.Length,data1$Sepal.Width,var.equal = T)
#Two sample t test
var.test(data1$Sepal.Length,data1$Sepal.Width)
#F test to compare two variances
var.test(data1$Sepal.Length,data1$Sepal.Width,alternative = "greater")
var.test(data1$Sepal.Length,data1$Sepal.Width,alternative = "less")
x<-seq(-3,3,.1)
y<-dnorm(x)
plot(x,y,type = "l")
plot(x,y)
qnorm(.95)
#normal table
qnorm(.95, lower.tail = F)
#normal table right tail
pnorm(1.96,lower.tail = F)
#probabilty
mad(data1$Sepal.Length,center = median(data1$Sepal.Length))
#Mean absolute deviation by default median theke deviation
mad(data1$Sepal.Length,center = mean(data1$Sepal.Length))
mad(data1$Sepal.Length)
library(psyche)
library("MASS", lib.loc="C:/Program Files/R/R-3.3.2/library")
describe(data1$Sepal.Length)
describe.by(data1$Sepal.Length)
describe.by(data1$Sepal.Length,group = data1$Species)
describe.by(data1$Sepal.Length,group = data1$Sepal.Length<=5.0)
describeBy(data1$Sepal.Length,data1$Species)
print(xy<-describe.by(data1$Sepal.Length,group = data1$Species))
xy$setosa$mad
data2<-read.table(file = "g:/fungal.txt",header = T)
data2
cor(data2$age,data2$clears)
cor.test(data2$age,data2$clears)
data2$severity<-factor(data2$severity)
attach(data2)
print(t<-table(treat,sex))
chisq.test(t)
fisher.test(t)
#odds ratio=(success/failure),1 means success=failure, 1< means success is more than failure
xtabs(~treat+sex+age)
age
detach(data2)
mydata<-read.table("location", header = T, sep = ",",row.names = "idofpatient")



#_____________Thesis____________

library(ggplot2)
library(xlsx)
library(openxlsx)
cl<-read.csv(file.choose(),header = T)
write.csv(dtc,"E:/222/CSE490/clinvar-conflicting/cli.csv")
library(foreign)
write.xlsx(dtc,"E:/222/CSE490/clinvar-conflicting/cli.xlsx")

tree_func <- function(final_model, 
                      tree_num) {
  
  tree <- randomForest::getTree(final_model, 
                                k = tree_num, 
                                labelVar = TRUE) %>%
    tibble::rownames_to_column() %>%

    mutate(`split point` = ifelse(is.na(prediction), `split point`, NA))
 
  graph_frame <- data.frame(from = rep(tree$rowname, 2),
                            to = c(tree$`left daughter`, tree$`right daughter`))
  
  graph <- graph_from_data_frame(graph_frame) %>%
    delete_vertices("0")

  V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
  V(graph)$leaf_label <- as.character(tree$prediction)
  V(graph)$split <- as.character(round(tree$`split point`, digits = 2))
  
  plot <- ggraph(graph, 'dendrogram') + 
    theme_bw() +
    geom_edge_link() +
    geom_node_point() +
    geom_node_text(aes(label = node_label), na.rm = TRUE, repel = TRUE) +
    geom_node_label(aes(label = split), vjust = 2.5, na.rm = TRUE, fill = "white") +
    geom_node_label(aes(label = leaf_label, fill = leaf_label), na.rm = TRUE, 
                    repel = TRUE, colour = "white", fontface = "bold", show.legend = FALSE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18))
  
  print(plot)
}

table(cll$CLASS)
1053800/32203
dim(dtc)

l()
rm(list = ls())

View(cl)
colnames(cl[,c(1,2,3,4,8,10,14,16,17,19,21,22,23,24,26,29,30,31,32,37,38,43,44,45)])
dtc<-cl[,c(1,2,3,4,8,10,14,16,17,19,21,22,23,24,26,29,30,31,32,37,38,43,44,45)]
View(dtc)
table(dtc$REF,dtc$ALT)
as.factor(dtc$ALT)
qplot(ALT, REF, data = dtc)

as.factor(dtc$CLNDN[dtc$CLNDN!="Deletion"])
as.factor(dtc$CLNDN)
x<-c(1,2,3)


fCLNDN<-na.omit((dtc$CLNVC))
fCLNDN2<-na.omit((dtc$CLNVC[dtc$CLNVC!="Deletion"]))

length(fCLNDN)

length(fCLNDN2)

cli<-read.xlsx("E:/222/CSE490/clinvar-conflicting/Final_Dataset.xlsx")
View(cli)
l1<-lm(cli$CLASS~cli$LoFtool+cli$CADD_PHRED+cli$CADD_RAW)
l2<-glm(cli$CLASS~cli$LoFtool+cli$CADD_PHRED+cli$CADD_RAW, family= "binomial")
fit.glm<-glm(CLASS~CHROM+CLNVC+IMPACT+cDNA_position+CDS_position+Protein_position+LoFtool+CADD_PHRED+CADD_RAW
             ,data=cli,family = "binomial")
class(cl$CHROM)


fit.glm<-glm(CLASS~CLNVC+IMPACT+cDNA_position+CDS_position+Protein_position+LoFtool+CADD_PHRED+CADD_RAW
             ,data=cll,family = "binomial")
plot(fit.glm$residuals)
fit.glm$fitted.values
fit.glm$effects
fit.g1<-glm(CLASS~cDNA_position+CDS_position+Protein_position+LoFtool, data= cli, family = "binomial")
c<-exp(fit.glm$coefficients)
c<-exp(fit.g1$coefficients)
fit<-glm(CLASS~CLNVC+IMPACT+cDNA_position+CDS_position+Protein_position+LoFtool+CADD_PHRED+CADD_RAW
             ,data=train,family = "binomial")
p<-exp(p)
p[p<.5]<-0
p[p>.5]<-1

p <- predict(fit.glm, test)
p<-exp(p)
p[p<.5]<-0
p[p>.5]<-1
summary(fit.glm)
View(p)
table(cl$CLASS)
table(p,test$CLASS)

library(openxlsx)
write.xlsx(cl,"e:/cl.xlsx")
install.packages("caTools")

library(caTools)
cll<-cl[,c(7,10,12,14,15,16,18,19,20)]
colnames(cl)
cll$CLASS<-as.factor(cll$CLASS)
cll$IMPACT<-as.factor(cll$IMPACT)
cll$CLNVC<-as.factor(cll$CLNVC)
class(cll$CLNVC)



set.seed(45)
split = sample.split(cl$CLASS, SplitRatio= 1/5)
training_set= subset(cl, split==T)
test_set=subset(cl, split==F)

rm(test_set)

set.seed(123)
ind<-sample(2,nrow(cll), replace = T, prob = c(.8,.1))
train<-cll[ind==1,]
test<-cll[ind==2,]


set.seed(45)

split = sample.split(cli$CLASS, SplitRatio = 1/2)
training_set= subset(cli, split==T)
test_set=subset(cli, split==F)

set.seed(45)

split = sample.split(cli$CLASS, SplitRatio= 1/2)
training_set= subset(cli, split==T)
test_set=subset(cli, split==F)

regressor<- lm(formula = hp~mpg, data = training_set)
p <- predict(regressor, test_set)
View(p)

library(randomForest)

library(dplyr)
library(ggraph)
library(igraph)
install.packages("caret")

library(caret)
install.packages("rlang")
library(rlang)

myLittleForest <- randomForest(CLASS~as.factor(CLNVC)+as.factor(IMPACT)+cDNA_position+CDS_position+Protein_position+LoFtool+CADD_PHRED+CADD_RAW,data=training_set)

myLittleForest <- randomForest(as.factor(train$CLASS)~as.factor(train$CLNVC)+as.factor(train$IMPACT)+train$cDNA_position+train$CDS_position+train$Protein_position+train$LoFtool+train$CADD_PHRED+train$CADD_RAW)
print(myLittleForest)

myLittleForest <- randomForest(CLASS~.,data=train)
print(myLittleForest)


myLittleForest <- randomForest(CLASS~cDNA_position+CDS_position+Protein_position+LoFtool,data=training_set)

attributes(myLittleForest)

myLittleForest$confusion

p1<-predict(myLittleForest,train)
confusionMatrix(p1, as.factor(train$CLASS))

head(p1)
head(training_set$CLASS)
fix(test_set)


#p2<-predict(myLittleForest,test)

plot(myLittleForest) #100

t<-tuneRF(train[,-2],train$CLASS,
          stepFactor = 0.5,
          plot = T,
          ntreeTry = 100,
          trace = T,
          improve = 0.05)



test_set[complete.cases(test_set),]

confusionMatrix(p2, as.factor(test_set$CLASS))

head(test_set)
dim(test_set)
dim(training_set)

plot(myLittleForest)


tree_func(myLittleForest,1)

myLittleForest <- randomForest(CLASS~CADD_PHRED+CADD_RAW,data=cli)


#PCA
pca1<-princomp(cl[,c(10,12,14,15,16,18,19,20)],scores = T, cor = T)
summary(pca1)
pca1$loadings
write.xlsx(as.data.frame(pca1$loadings),"e:/pca1loadings.xlsx")
loadings(pca1)
plot(pca1)
screeplot(pca1,type = "line",main = "screeplot")
biplot(pca1)
pca1$scores[1:10,]
fa1<-factanal(cl[,c(10,12,14,15,16,18,19,20)],factors = 3)
fa1
fa1<-factanal(cl[,c(10,12,14,15,16,18,19,20)],factors = 3,rotation = "varimax")
fa1




# Print the summary of our model.
print(myLittleForest)
print(importance(myLittleForest, type=2))
plot(myLittleForest)

l2$control
library(haven)

cl<-read_sav("E:/222/CSE490/clinvar-conflicting/ClinVar.sav")
View(cl)


summary(l1)
summary(l2)
create_report(cli)



X~binomial(10,3)
dbinom(2,10,0.3)
#P(X=2)
pbinom(3,10,0.3)
#P(X<=3)
pbinom(2,10,0.3)
#P(X<3)
1-pbinom(4,10,0.3)
#P(X>4)
qbinom(0.5, 10, 0.3)  
# Binomial Example2 : X ~ binomial(10, 0.3)
# Find value k such that P(X <= k) = 0.5 

pnorm(z<-(16-20)/sqrt(4),lower.tail = F)
pnorm(16,20,sqrt(4),lower.tail = F)

pt(-2.15,19)
pt(-2.15,19,lower.tail = F)


pt(1.96,100,lower.tail =F)*2
pnorm(1.96,lower.tail = F)*2
pchisq(3.841,1,lower.tail = F)
pf(2.5,df1 = 2,df2 = 15,lower.tail = F)


x<-rnorm(1000,10.06207,sqrt(0.11677))

l<-(x<=9.5)
l1<-(x>=10.5)
n<-length(x)
lp1<-(sum(l1)/n)*100
lp<-sum(l)/n *100


bnx<-c(1:1)
bny<-c(1:1)
for(i in 1:1000)
{
  if(x[i]<=9.5){
    c(bnx[i]<-x[i])
  }
  if(x[i]>=10.5){
    c(bny[i]<-x[i])
  }
}

(length(bnx<-na.exclude(bnx))/1000)*100
(length(bny<-na.exclude(bny))/1000)*100

x11<-x[x>=10.5]
x22<-x[x<=9.5]
(length(x11)/1000)*100
(length(x22)/1000)*100


x

data("mtcars")
attach(mtcars)
ls()
head(mtcars,3)
mtcars
nrow(mtcars)
sample(1:32,10,replace = F)
#First avoid the missings and then sampling
ndat<-mtcars[sample(1:32,10,replace = F),]
data(iris)
boxplot(iris[,1])
boxplot(iris[,1]~iris[,5])
par(mfrow=c(2,2))
boxplot(iris[,2]~iris[,5])
for(i in 1:4)
{
  boxplot(iris[,i]~iris[,5])
}

bmp("g:/test.bmp")
par(mfrow=c(2,2))
for(i in 1:4)
{
  boxplot(iris[,i]~iris[,5])
}
dev.off()
#for graphical files
par(mfrow=c(1,1))
mpg
l<-lm(mpg~disp+hp)
l<-lm(mpg~-1+disp+hp)
summary(l)
vcov(l)
v<-diag(vcov(l))
sqrt(v)
coefficients(l)
confint(l,level = 0.95)
fitted(l)
residuals(l)
anova(l)
l$fitted
l$residuals
rss<-var(l$fitted.values)*31
#computer by deafault n-1 kore
ess<-var(l$residuals)*31
tss<-ess+rss
help(package=rattle)
rattle()
install.packages("https://cran.r-project.org/bin/windows/contrib/3.3/RGtk2_2.20.31.zip", repos=NULL)
library(Rcmdr)

m<- matrix(c(150,50,240,60),2,2)
print(t<-table(data2$X1,data2$X2))
chisq.test(m)
 banknote
library(rattle)
rattle()
install.packages("mclust")
banknote
write.table(banknote,file='d:/banknote.csv',sep=",")

m_factor<-factor(m_vector,ordered = T, levels = c("low","medium","high"))

cost_2014<-c(8.6,8.5,8.1)

cost_2014>8.3
cost_2014[cost_2014>8.3]

movie_vector<-c("Naimul", "Hasan","Helal","New")
movie_array<- array(movie_vector,dim=c(4,3))
movie<-list(name="Toy Story",
            year=1995,
            genre=c("Animation","Adventure","Comedy")
)
movies<-data.frame(name=c("Naimul", "Hasan","Helal","New"),
                   year=c(1001:1004))
movies['length']<-c(12,51,31,41)
movies<rbind(movies, c(name="md",year=45,length=96))
movies["length"]<-NULL

%in% #Found in
myFunction<-function(){
  y<<-3.14
  temp<-'Hello Nsimul'
  return(temp)

}
#Here y is an global variable
age<-c(120,15,61,45)
integer_age<-as.integer(age)

tryCatch(
  for(i in 1:3)
  {
    print(i+'a')
  }
  , error = function(e){
    print("Found Error.")
  }
)

tryCatch(
  as.integer("A"),
  warning= function(e){
    print("Warning")
  }
)

read.csv("location")
install.packages("readxl")
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
plot(TF1$ForeignEx, TF1$Total,main= "Corelation Between ",xlab="Foreign Ex",ylab="Total Tourist", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
plot(TF2$ForeignEx, TF2$Total,main= "Corelation Between ",xlab="Foreign Ex",ylab="Total Tourist", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
plot(TFT$ForeignEx, TFT$Total,main= "Corelation Between ",xlab="Foreign Ex",ylab="Total Tourist", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
cor(TFT$ForeignEx, TFT$Total)
cor(TF1$ForeignEx, TF1$Total)
cor(TF2$ForeignEx, TF2$Total)
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
TF0<-TF1[2:15,]
t.test(TF0$ForeignEx,TF2$ForeignEx,var.equal=T,alternative = "less")

tf1<-lm(TF1$Total~TF1$Year)
summary(tf1)
#R-sq: 0.7865
cv(TF1$Total)*(sqrt(1-0.7865))

tf2<-lm(TF2$Total~TF2$Year)
summary(tf2)
#R-sq: 0.5429
cv(TF2$Total)*(sqrt(1-0.5429))



read_excel("location")
data()#to see all the datasets
text<-readLines("location")#by line
text<-scan("location", "")#by word
nchar(text)
file.size("location") #in bytes

hw<-c("Hello","Naimul")
paste(hw, collapse = " ")
#for concat strings

as.Date("2020-01-01")
library(stringr)
m<- matrix(c(1:6),2,3)
write(m, file = "location.txt", ncolumns = 3, sep = " ")
write.csv(df, file= "location.csv", row.names = F)
write.table(df, file = "location.csv",row.names = F, col.names = F, sep = ",")
dt <- read.spss("E:/Link_Final.sav")
dt <- as.data.frame(dt)
View(dt)
fix(dt)
td<-dt[,-13]
?str_split_regex()
length(td[,1])

write.csv(x = as.factor(unlist(dat)), file = "E:/f.csv")
unique(trimws(as.character(x),"r"))
x<-trimws(unlist(dat),"r")
x<-trimws(x,"l")
trimws(x)
x<-gsub("[.]","",x)



dat<-(str_split(td$Others,"\\+"))
View(dat)
fix(dat)
td$id<-c(1:length(td[,1]))

df<-data.frame(id = factor(rep(td$id,times= lengths(dat)),levels =td$id),
               StriN = x)
df1<-as.data.frame(cbind(id=td$id,
                         SubCategory=td$SubCategory,Rank=td$Rank,
                         CourseName=as.character(td$CourseName),Viewers=td$Viewers,
                         Student=td$Student,Professor=td$Professor,
                         Accountant=td$Accountant,Salesperson=td$Salesperson,
                         ITConsultant=td$ITConsultant,
                         Teacher=td$Teacher,
                         SoftwareDeveloper=td$SoftwareDeveloper,
                         table(df$id,df$StriN)))
write.csv(df2, file = "E:/f.csv")
write.foreign(df2, "E:/f.txt", "E:/f.sps",   package="SPSS")
install.packages("xlsx")

library(xlsx)
write.xlsx(df2,"E:/finalLink.xlsx")

df2<-data.frame(td[,c(1:11)],df1[,c(15:122)])
View(df2)

df2<-as.data.frame(table(df$id,df$StriN))
x<-"  .Some text. "
trimws(x)
gsub("[.]","",x)


View(df)
View(td)
View(df1)
fix(df1)
as.factor(colnames(df1))
t1<-table(df$id,df$StriN)
t2<-table(df$id,df$StriN)
View(t1)
View(t2)

write.xlsx(df, file = "location.xlsx", sheetName = "Sheet1", col.names = T, row.names = F)
naim <- read.csv("C:/Users/Dolphin-Pc/Desktop/mah.csv", sep =',')
View(naim)
dat<-str_split(naim$StriN,",")
df<-data.frame(id = factor(rep(naim$id, lengths(dat)),levels = naim$id),
               StriN = unlist(dat))
df1<-as.data.frame(cbind(id=naim$id,table(df$id,df$StriN)))
library(car)
df1$Bhai<-recode(df1$Bhai,"1=40")
df1$Naimul<-recode(df1$Naimul,"1=20")


save()

toupper()
tolower()
chartr("a","_",summary1[1])
#here "a" will be replaced by "_" from first summary1 line

char_list<- strsplit(summary1[1]," ")
word_list<- unlist(char_list)

sub_string<- substr(summary1[1], start = 4, stop = 50)
trimws(sub_string)

library(stringr)
str_sub(summary1[1],-8,-1)

actors.birthday<-as.POSIXct(bestActors$Date.of.Birth, origin="1970-01-01")
actors.birthday<-as.Date(actors.birthday)
#From all seconds to the date format conversion
actors.birthday<- as.Date(bestActors$Date.of.Birth,"%Y/%m/%d")
#from factors to date conversion

as.Date("1999/06/27")-as.Date("1959/01/01")
#date operations
as.Date("1999/06/27")-14
Sys.Date()# can do seq or other operations
Sys.time()
date()
grep("@.*",c("afafa@eresf.com","sefaf@asfaf.com","fea@afaf.com"), value = T)
grep("@.+",c("afafa@eresf.com","sefaf@asfaf.com","fea@afaf.com","ad@"), value = T)
grep(".+@.*",c("afafa@eresf.com","sefaf@asfaf.com","fea@afaf.com"), value = T)

gsub("@.*","@yahoo.com",c("afafa@eresf.com","sefaf@asfaf.com","fea@afaf.com"))

matches<-regexpr("@.*\\.",email_df[,'Email'])
email_df[,'Domain']=regmatches(email_df[,'Email'], matches)

install.packages("SparkR")

saveRDS(df1, file = "location.Rda")
df2<-readRDS("location.Rda")
#SIngle Object File

save(df1, df2, file="multi.Rda")
OR
save(list=c("df1","df2"), file="multi.Rda") 
load("multi.Rda")
#Multiple Object File

save.image("all.Rda")
load("all.Rda")
# for all data structures in the workspace


library(foreign)
dtx<-read.dta(file = "d:/BDHS14.dta",convert.factors = T)
View(dtx)
adtx<-attributes(dtx)
adtx$var.labels
library(DMwR)
install.packages("DMwR")
dtx2<-dtx[,c(1,5)]
View(dtx2)
length(dtx2[dtx2$bmi>=9998,2])
length(dtx2[,2])
dtx2[dtx2$bmi>=5000 & dtx2$bmi<=6000,2]
dtx2[dtx2$bmi>=4000 & dtx2$bmi<=5000,2]
dvsn<-summary(dtx$division)
dtx[dtx$bmi>=4000 & dtx$bmi<=6000 & dtx$division=="Sylhet",5]
dtx1<-dtx[dtx$bmi<9998,]
View(dtx1)
dvsn1<-summary(dtx1$division)
sampledtx1<-dtx1[sample(1:nrow(dtx1),nrow(dtx1)*.2,replace=F),]
summary(sampledtx1$division)
rm(sampledtx1)
View(dtx2)
dtx2<-dtx1[order(dtx1$division),]
dvsn<-levels(factor(dtx2$division))

sp <- split(dtx2, list(factor(dtx2$division)))
samples <- lapply(sp, function(x) x[sample(1:nrow(x), nrow(x)*.20, FALSE),])
sdtx2 <- do.call(rbind, samples)
View(sdtx2)
summary(sdtx2$division)
#Stratified Sampling
summary(dtx2$division)


for(i in 1:7)
{
  c(s_d[i]<-sd(dtx2[dtx2$division==dvsn[i],5])) 
  c(mn[i]<-mean(dtx2[dtx2$division==dvsn[i],5]))
}
mn<-c(10,12)
s_d<-c(14,16)
mn
s_d
install.packages('car')
library(car)
dtx2$bmicat<-recode(dtx2$bmi,"1000:1849='underweight';1850:2299='Normal'; else='Overweight'")
tdtx2<-table(dtx2$bmicat,dtx2$plor)
chisq.test(tdtx2)
barplot(tdtx2,xlab="Urban or Rural",ylab = "Body Mass Index",
        col=c("green","blue","yellow"),
        legend.text = T,
        args.legend = list(x=8,y=1600))
tdtx2
bd<-dtx2
attach(bd)
boxplot(age, data = bd, main= "Box plot of age",xlab="From data BDHS in 2014",ylab="Age", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
count<-table(religion)
barplot(count, main= "Bar plot of religion",xlab="Several Religion",ylab="Number of people", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
box()
x<-table(division)
pie(x, data = bd, main= "Pie diagram of several division",cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
box()
hist(bmi,main= "Histogram of bmi",xlab=" Body Mass Index",ylab="Frequency", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
box()
boxplot(age~division, data = bd, main= "Box plot of several division",xlab="Divisions",ylab="Age", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)
box()

cor(age, bmi)

plot(age, bmi,main= "Corelation Between ",xlab="Divisions",ylab="Age", cex=1, cex.main=1.5, cex.lab=2.0, cex.axis=1.5, col=3, col.main=2, col.lab=4, col.axis=9)

mod <- lm(bmi~age)
summary(mod)
plot(mod)
write.csv(dtx,"d:/bdhs.csv")
write.table(dtx,"d:/bdhs.csv",sep=",")
install.packages("Hmisc")
library(Hmisc)
install.packages("setRNG")
library("setRNG")
rng(0,'twister')
a=5
b=500
y=a*randn(1000,1)+b

rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }
r <- rnorm2(100,0,2)
mean(r)
sd(r)   
x[x>=10]
sum(x[x>=10])

s<-12;u<-5;s%%u#Modoulus
s<-20;u<-3;s%/%u#Integer Division

r2<-rnorm(100,0,2)
print(meansd(r2))

r1<-2*scale(rnorm(100))

y<-cut(x2,5)
y
y<-x2
y

banknote
vif(banknote)
help(vif)
library("Hmisc")
install.packages("fmsb")
library("fmsb")


for(i in 1:10)
{
  x<-rnorm(25)
  print(t.test(x))
}
meannn<-function(x,n)
{
  sum((1/(n-1))*(x-mean(x))^2)
}
x<-c(1:10)
meannn(x,10)
____________________________RJDBC________________________________
install.packages("RJDBC")
library(RJDBC)
drv <- JDBC("com.ibm.db2.jcc.DBZDriver", "db2jcct4.jar", identifier.quote="'") 

?anova()
X<-matrix(c(750,250,350,650),2,2,byrow = T)

ch<-chisq.test(X)
ch
n<-100000
x<-seq(from=0, to=1, by=.0001)
x
x
for(i in 1:10000){
  (c(nai[i]<-exp(exp(x[i]))))
}
nai
sum(nai)
set.seed(1234)
u<-runif(100000)
v<-exp(exp(u))
tv<-mean(v)
tv
se<-sqrt(var(v)/n)
se
set.seed(1200)
u<-runif(100000)
v<-exp(exp(u))
tv<-mean(v)
tv
se<-sqrt(var(v)/n)
se

set.seed(1200)
u<-runif(100000)
v<-(1-u^2)^(3/2)
tv<-mean(v)
tv
se<-sqrt(var(v)/n)
se
(3*3.1416)/16

set.seed(1234)
u<-runif(100000)
v<-(1-u^2)^(3/2)
tv<-mean(v)
tv
se<-sqrt(var(v)/n)
se

set.seed(1234)
u<-runif(100000,-2,2)
v<-exp(x+x^2)
tv<-mean(v)
tv
se<-sqrt(var(v)/n)
se
x1<-(1/u^2)
x2<-((1/u)-1)
x3<-((1/u)-1)^2

set.seed(1234)
u<-runif(100000)
v<-x1*x2*(x3+1)^(-2)
tv<-mean(v)
tv
se<-sqrt(var(v)/n)
se
xx<-exp(-((1/u)-1)^2)
x1<-1/(u^2)

set.seed(1234)
u<-runif(100000)
v<-x1*xx
tv<-mean(v)
tv
se<-sqrt(var(v)/n)

set.seed(1234)
u<-runif(100000)
x<-runif(100000)

v<-exp((x+u)^2)
tv<-mean(v)
tv
se<-sqrt(var(v)/n)
se

set.seed(1234)
u<-runif(100000)
v<-(u-.5)*(exp(u))
tv<-mean(v)
tv
se<-sqrt(var(v)/n)
se

set.seed(90)
t<-runif(min = 0, max= 1, n=10000)
mean(exp(exp(t)))

set.seed(100)
u<-runif(10000)
u1<-(-2)+4*u
mean(exp(u1+u1^2))
se<-sqrt(var(exp(u1+u1^2))/10000)
se
#question 1
r1<-rbinom(n= 100, size= 1,prob = 1/3)
l<-length(r1[r1==1])
p<-(l/100)
p

r2<-rbinom(n= 1000, size= 1,prob = 1/3)
l<-length(r2[r2==1])
p2<-(l/1000)
p2

r3<-rbinom(n= 10000, size= 1,prob = 1/3)
l<-length(r3[r3==1])
p3<-(l/10000)
p3
#question_1.2
r1<-rbinom(n= 100, size= 1,prob = 1/4)
l<-length(r1[r1==1])
p<-(l/100)
p

r2<-rbinom(n= 1000, size= 1,prob = 1/4)
l<-length(r2[r2==1])
p2<-(l/1000)
p2

r3<-rbinom(n= 10000, size= 1,prob = 1/4)
l<-length(r3[r3==1])
p3<-(l/10000)
p3

#ctrl+ Alt+ Click
#quest_3
u<-runif(10000)
for(i in 1:10000)
{
  if (u[i]<0.35) 
  {
    c(x[i]<-3)
  }
  else if(u[i]<0.65)
  {
    c(x[i]<-1)
  }
  else if(u[i]<0.85)
  {
    c(x[i]<-2)
  }
  else
  {
    c(x[i]<-4)
  }
}
l<-length(x[x==1])
p<-l/10000
p
l<-length(x[x==2])
p<-l/10000
p
l<-length(x[x==3])
p<-l/10000
p
l<-length(x[x==4])
p<-l/10000
p
x<-c(1:9)
for(i in 1:99)
{
  c(x[i]<-1/(i*(i+1)))
}
sum(x)
rm(x)

u<-runif(10000)
t<-(1/u-1)/(2*(u^2)-2*u+1)
mean(t)

u<-runif(100000)
t<-exp(-((1/u)-1)^2)/u^2
2*mean(t)

x<-runif(100000)
y<-runif(100000)
f<-exp((x+y)^2)
mean(f)

for(i in 1:99)
{k[i] <- (1/((i)*(i+1)))}
sum(k)

library(leaflet)
library(htmlwidgets)
library(IRdisplay)
lower_lon = 91.847
upper_lon = 91.849
lower_lat = 24.90
upper_lat = 24.92
center_lon = (lower_lon + upper_lon)/2
center_lat = (lower_lat + upper_lat)/2
subset <- sh 
road_map <- leaflet(subset) %>% 
  setView(center_lon,center_lat, 17)%>% 
  addProviderTiles(providers$Thunderforest.TransportDark)%>%
  addMiniMap(
    tiles = providers$Esri.WorldImagery,
    toggleDisplay = TRUE)%>% 
  addCircleMarkers(lng = subset$Longitude, 
                   lat = subset$Latitude, 
                   #popup = subset$Stn_Name,
                   fillColor = "Black", 
                   fillOpacity = 1, 
                   radius = 4, 
                   stroke = F) 
road_map
saveWidget(road_map, file="G:/ML0151/road_map1.html", selfcontained = F) 

install.packages("survival")
install.packages("ggplot2")
options(repos='http://cran.rstudio.com/')
library(survival)
library(ggplot2)
lg<-lung
head(lg)
View(lg)
recodest<-function(x){
  if(x==1){rs=0}
  if(x==2){rs=1}
  return(rs)
}
for(i in 1:length(lg$status)){
  lg$rcStatus[i]<-recodest(lg$status[i])
}
mySurv<-Surv(lg$time,lg$rcStatus)
head(mySurv)
myfit<-survfit(mySurv~lg$sex)
myfit
summary(myfit)

plot(myfit)
plot(myfit, conf.int = "none")
abline(h=0.5)
abline(v=median(lg$time))
abline(v=310)
mf1<-survfit(mySurv~lg$sex)
plot(mf1)
table(lg$sex)
plot(mf1, col = c("red","blue"),mark=3)
plot(mf1,conf.int = "none", col = c("red","blue"),xlab = "time",ylab = "Cumulative Survival Probabilty")

legend("topright",c("male","female"),col=c("red","blue"),lty = 1)
abline(h=0.5)
abline(v=270,col= "blue")
abline(v=426,col= "red")
?survdiff()
survdiff(mySurv~lg$sex)
plot(mf1,fun = "event",conf.int = "both", col = c("red","blue"))
summary(mf1)
survfit(Surv(lg$time,lg$status)~lg$sex)
View(lg)

library(foreign)
read.csv("G:/STA436/lab2.csv",sep=",")
sl<-read.csv("G:/STA436/lab2.csv",sep=",")
fix(sl)
head(sl)

myS<-Surv(time = sl$time, event = sl$status)
head(mySurv)
m<-survfit(myS~sl$group)
mt<-survdiff(myS~sl$group, rho = 0)
mt<-survdiff(myS~sl$group, rho = 5)
?survdiff

m
plot(m,conf.int = "none", col = c("red","blue"),xlab = "time",ylab = "Cumulative Survival Probabilty", mark = 3)

legend("topright",c("sample-1","sample-2"),col=c("red","blue"),lty = 1)
summary(m)

plot(myfit, conf.int = "none")
abline(h=0.5)
abline(v=median(lg$time))
abline(v=310)
mf1<-survfit(mySurv~lg$sex)
plot(mf1)
table(lg$sex)
plot(mf1, col = c("red","blue"),mark=3)
plot(mf1,conf.int = "none", col = c("red","blue"),xlab = "time",ylab = "Cumulative Survival Probabilty")

legend("topright",c("male","female"),col=c("red","blue"),lty = 1)
abline(h=0.5)
abline(v=270,col= "blue")
abline(v=426,col= "red")
?survdiff()
survdiff(mySurv~lg$sex)
plot(mf1,fun = "event",conf.int = "both", col = c("red","blue"))
summary(mf1)
survfit(Surv(lg$time,lg$status)~lg$sex)
View(lg)

bio<-read.csv("E:/222/STA436/project_436_2018.csv",sep = ",", header = T)
View(bio)
fix(bio)
mod <- lm(bio$bmi ~ bio$AGE+bio$WAIST+bio$ARM)
colnames(bio)
?sapply()
install.packages("tables")
#EDA
library(tables)
install.packages("DataExplorer") 
library(DataExplorer)
?DataExplorer

create_report(bio)


bio$bmicat[bio$bmicat==1]<-"Underweight"
bio$bmicat[bio$bmicat==2]<-"Normal"
bio$bmicat[bio$bmicat==3]<-"Overweight"
bio$bmicat[bio$bmicat=="Obese"]<-NA
library(plyr)
bio<-bio[-119,]
bio<-bio[-16,]
count(bio$bmicat=="Obese")
bio[bio$bmicat=="Obese",]

ggplot(data = bio, mapping = aes(x = bio$bmicat, y = bio$bmi)) +
  geom_boxplot()
#1
ggplot(data = bio, mapping = aes(x = bio$bmicat, y = bio$AGE)) +
  geom_boxplot()
#2
bio$fCLASS<-as.factor(bio$CLASS)

ggplot(data = bio, mapping = aes(x = bio$fCLASS, y = bio$AGE)) +
  geom_boxplot()
#3

ggplot(data = bio, mapping = aes(x = bio$fCLASS, y = bio$bmi)) +
  geom_boxplot()
#4
install.packages("hexbin")
library(hexbin)
ggplot(data = bio) +
  geom_hex(mapping = aes(x = bio$fCLASS, y = bio$bmi))

#5
barp<-ggplot(bio,aes(x=1,y=sort(bio$bmicat),fill=sort(bio$bmicat))) +
  geom_bar(stat = "identity")
print(barp)
print(barp<-barp+coord_polar(theta = 'y'))
print(barp<-barp+theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.y=element_blank(),
  panel.background = element_blank())+
    labs(y="BMI Catagory")
)
#7
View(bio)
library(ggplot2)
install.packages("GGally")
library(GGally)
par(mfrow=c(5,1))

hist(bio$bmi[bio$fCLASS=="1"],main = "Histogram For Class 1")
hist(bio$bmi[bio$fCLASS=="2"],main = "Histogram For Class 2")
hist(bio$bmi[bio$fCLASS=="3"],main = "Histogram For Class 3")
hist(bio$bmi[bio$fCLASS=="4"],main = "Histogram For Class 4")
hist(bio$bmi[bio$fCLASS=="5"],main = "Histogram For Class 5")
dev.off()
install.packages("easyGgplot2")
install.packages("devtools")

library(devtools)
install_github("kassambara/easyGgplot2")
library(easyGgplot2)
ggplot2.histogram(data=bio, xName='bmi',
                  groupName='fCLASS', legendPosition="top")

ggplot2.histogram(data=bio, xName='bmi',
                  groupName='fCLASS', legendPosition="top",
                  alpha=0.5, addDensity=TRUE,
                  addMeanLine=TRUE, meanLineSize=.5)


ggplot(bio,aes(x=bio$AGE,y=bio$bmi,color = bio$bmicat)) + geom_point(shape=19) + xlab("Age of Students") + ylab("Body Mass Index") +
  labs(colour = "BMI Catagory") + ggtitle("Scatterplot")
data(flea)
data(flea)
View(flea$species)
View(bio)

View(bio$bmicat)

bio$bmicat<-as.factor(bio$bmicat)

ggscatmat(flea, columns = 2:4, color="species", alpha=0.8)
ggscatmat(bio[,c(4,6,7,8,9,11)], columns = 1:6), color = "bio$bmicat",alpha = 0.8)
pairs(bio[,c(4,6,7,8,9,11)], col=bio$bmicat)



?ggscatmat

?create_report
create_report(bio)

summary(bio$SCHOOL)
table_freq(bio$SCHOOL)
plot_correlation(bio[,c(4,6,7,8,9,11)])

ggplot(data = bio, mapping = aes(x = bio$GENDER)) +
  geom_histogram(binwidth = 0.01)
q111<-qplot(bio$bmicat, 
            geom = "bar",
            fill= I("Pink"),
            colour=I("Blue"),
            xlab = "BMI Cata",
            ylab= "Number of Vehicles",
            main= "Cylinders in mtcars"
)


hist(1/bio$bmi^2)
shapiro.test(1/bio$bmi^2)
bio$bmit<-1/bio$bmi^2
f<-lm(bio$bmit~bio$AGE)
summary(f)
sqrt(1/2.603e-04)
f1<-lm(bio$bmi~bio$AGE)
summary(f1)
summary()
summary(mod)
shapiro.test(bio$bmi)
shapiro.test(bio$AGE)
shapiro.test(bio$WAIST)
?jtrans()
install.packages("jtrans")
library(jtrans)
shapiro.test(jtrans(bio$bmi))
shapiro.test(mod$residuals)

mean(mod$residuals)
qqnorm(mod$residuals)
install.packages("HH")

?shapiro.test


x<-c(10,15,16,18,20,25,18,15,16,21)
mean(x)
sd(x)

se<-sd(x)/10^.5
set.seed(4589)
a<-replicate(100, mean(sample(x,10,replace =T)))
mean(a)
sd(a)
hist(a)

set.seed(4589)
a<-replicate(10000, mean(sample(x,100,replace =T)))
mean(a)
sd(a)
hist(a)
print(se<-sd(a)/100^.5)

shapiro.test(a)

length(a)

mean(a)+ qt(0.975,100-1)*sd(a)*c(-1,+1)/sqrt(100)

m<-mean(a)

se<-sd(a)/100^.5

a<-replicate(10000, mean(sample(a,100,replace =T)))

c(m,mean(a))
m

quantile(a, probs = c(0.025,0.975))

#________Sampling Technique II Lab____________________
#5.1
n<-c(1:10)
x<-c(50,30,45,25,40,26,44,35,28,27)
dt<-data.frame(n,x)
dt$cum<-cumsum(x)
z<-c(1:3)

max(dt$cum)
head(dt)
b<-c(1:3)
for(i in 1:10)
{
  #set.seed(45654)
  y<-round(runif(10,min = 1, max = max(dt$cum)),0)
  #print(y)
  for(j in 1:10)
  {
    if(y[i]<dt$cum[j])
    {
      z[i]<-dt$n[j]
      break
      
    }
  }
  b<-unique(z)
  #b<-(sample(unique(z),4, replace = F))
  
}
z
b[1:4]

z
y
cumsum(1:10)
#cummin(c(3:1, 2:0, 4:2))
CTM<-function(u,m,k,l)
{  
  f<-cumsum(u)
  o<-length(m)
  for(i in 1:o)
  {
    set.seed(k)
    y<-round(runif(o,min = 1, max = max(f)),0)
    for(j in 1:o)
    {
      if(y[i]<f[j])
      {
        z[i]<-m[j]
        break
        
      }
    }
    b<-unique(z)
  }
print(b[1:l])
}
CTM(dt$x,dt$n,2013134027,4)
#CTM(Cluster Size,cluster number, seed, sample size)


for(i in 1:length(dt$n))
{
  y<-round(runif(1,min = 1, max = length(dt$n)),0)
  k<-round(runif(1,min = 1, max = max(dt$x)),0)
  if(k<=dt$x[y])
  {
    z[i]<-y
  }
}

LM<-function(g,h,u,l)
{
  set.seed(u)
  for(i in 1:length(h))
  {
    y<-round(runif(1,min = 1, max = length(h)),0)
    k<-round(runif(1,min = 1, max = max(g)),0)
    if(k<=dt$x[y])
    {
      z[i]<-y
    }
  }
  b<-unique(z)
  print(b[1:l])
}  

LM(dt$x,dt$n,345,4)
#LM(Cluster Size, Cluster Number,Seed, Sample Size)
#9A.2 Cochran
n<-c(1:7)
x<-c(3,1,11,6,4,2,3)
dt1<-data.frame(n,x)
LM(dt$x,dt$n,2013134027,2)


#5.3

area<-c(5.2,5.9,3.9,4.2,4.7,4.8,4.9,6.8,4.7,5.7,5.2,5.2,4.9,4.0,1.3,7.4,7.4,4.8,6.2,6.2)
yield<-c(28,29,30,22,24,25,28,37,26,32,25,38,31,16,6,61,61,29,47,47)
cdt<-data.frame(area,yield)
sum(cdt$area)

p<-area/484.5

(484.5/(20*100))*sum(yield/area)
sum((yield/area)^2)
sum(p)


yps<-(1/(20*100))*sum(yield/(area/484.5))
sum(area/484.5)


vpps<-(1/(20*19))*(((484.5^2)*sum((yield/area)^2))/(100*100)-20*(yps^2))
sqrt(vpps)
vsr<-(100/20)*{((sum((yield^2)/(area/484.5)))/(20*100*100))-(1/100)*((yps^2)-vpps)}

v_pps<-function(x,y,N,n,tx)
{
  yps<-(1/(n*N))*sum(y/(x/tx))
  vpps<-(1/(n*(n-1)))*(((tx^2)*sum((y/x)^2))/(N^2)-n*(yps^2))
  svpps<-sqrt(vpps)
  vsr<-(N/n)*{((sum((y^2)/(x/tx)))/(n*N*N))-(1/N)*((yps^2)-vpps)}
  prec<-((vsr-vpps)/vpps)*100
  
  return(list("Y_bar_PPS"=yps,"Variance_PPS"=vpps,"Standard_Error_PPS"=svpps,"Variance_SRS"=vsr,"Precision_PPS"=prec))
}
v_pps(area,yield,100,20,484.5)
v_pps(al,yl,100,20,484.5)

sum((yield^2)/(area/484.5))
484.5*sum((yield^2)/area)
((vsr-vpps)/vpps)*100

al<-c(4.8,4.1,1.3,5.2,6.9,6.0,2.0,6.3,5.2,4.2,4.8,5.9,5.8,5.8,5.1,4.7,5.6,5.2,4.0,4.6)
length(al)
yl<-c(22,19,6,25,54,43,4,40,28,29,22,39,39,44,30,27,34,31,18,31)
length(yl)

q<-(1/20^2)*((sum((yield^2)/(area/484.5)*(100-(484.5/area)))))
q/vpps
ldt<-data.frame(al,yl)

#5.5

#a. Murthy Method
s<-c(1:8)
x<-c(50,30,25,40,26,44,20,35)
y<-c(60,35,30,44,30,50,22,40)

dt<-data.frame(s,x,y)
dt$p<-x/sum(x)
head(dt)
x1<-dt$x[5]
x2<-dt$x[7]
y1<-dt$y[5]
y2<-dt$y[7]
p1<-dt$p[5]
p2<-dt$p[7]

ym<-(1/(2-p1-p2))*(((y1/p1)*(1-p2))+((y2/p2)*(1-p1)))

vym<-(((1-p1)*(1-p2)*(1-p1-p2))/(2-p1-p2)^2)*((y1/p1)-(y2/p2))^2
svym<-sqrt(vym)

#b. Horvitz-Thompson

s<-sum(dt$p/(1-dt$p))

pp1<-p1*(s+1-(p1/(1-p1)))
pp2<-p2*(s+1-(p2/(1-p2)))
pp12<-p1*p2*((1/(1-p1))+(1/(1-p2)))

yht<-((y1/pp1)+(y2/pp2))
vyht<-(((pp1*pp2)-pp12)/pp12)*((y1/pp1)-(y2/pp2))^2
svyht<-sqrt(vyht)

# p1<-c(4.32,4.16,3.06,4,4.12,4.08,5.16,4.40,4.20,4.28)
# p2<-c(4.84,4.36,4.24,4.84,4.68,3.96,4.24,4.72,4.66,4.36)
# p3<-c(3.96,3.5,4.76,4.32,3.46,3.42,4.96,4.04,3.64,3)
# p4<-c(4.04,5,3.12,3.72,4.02,3.08,3.84,3.98,5,3.52)
# f<-c(1,2,3,4,5,6,7,8,9,10)
# wht<-data.frame(f,p1,p2,p3,p4)
# wht$t2<-sum(wht[,c(2)])
# View(wht)
# delete(wht$t2)
# 
# wht$t2<-p1+p2+p3+p4

#9.2

m1<-c(266,890,311,46,174,31,17,186,224,31,102,46,31,109,275,128,125,267,153,152,84,21,52,10,0,48,94,123,87,89,109,0,310,3)
m2<-c(129,57,64,11,163,77,278,50,26,127,252,194,350,0,572,149,275,114,387,53,34,150,224,185,157,244,466,203,354,816,242,140,66,590,747,147)
m3<-c(247,622,225,278,181,132,659,403,281,236,595,265,431,190,348,232,88,1165,831,120,987,938,197,614,187,896,330,485,60,60,1051,651,552,968,987)
m4<-c(347,362,34,11,133,36,34,61,249,170,112,42,161,75,68,0,247,186,473,0,143,198,65,0,308,122,345,0,223,302,219,120,199,35,0,0)
mean(102,105,200,88)
sum(m1)


x1<-(m1-mean(m1))
s1<-(1/33)*sum(x1^2)

x2<-(m2-mean(m2))
s2<-(1/35)*sum(x2^2)

x3<-(m3-mean(m3))
s3<-(1/34)*sum(x3^2)

x4<-(m4-mean(m4))
s4<-(1/35)*sum(x4^2)
si<-c(s1,s2,s3,s4)

#First Estimate
yEs<-round(((102*mean(m1))+(105*mean(m2))+(200*mean(m3))+(88*mean(m4)))/(4*124),0)

Mi<-c(102,105,200,88)
M<-Mi/124
y<-c(135,225,471,141)
x<-(M*y-290)^2

sb<-(1/3)*sum(x)
mi<-c(34,36,35,36)
vy1<-(M^2)*((1/mi)-(1/Mi))*si

vy<-(1/(4*12))*sum(vy1)

vybar<-(sb/6)+vy
sevybar<-sqrt(vybar)

#Second Estimate
yibar<-mean(y)
yid<-(y-yibar)^2
nsb<-(1/3)*sum(yid)
vyi<-((1/mi)-(1/Mi))*si

vyibar<-(((1/4)-(1/12))*nsb)+(1/48)*sum(vyi)
svyibar<-sqrt(vyibar)

#Quantile / Percentile
quantile(dt$bmi,c(.30,.50,.75))

#____________Yousuf Start_____
library(foreign)
dt <- read.spss("E:/222/STA436/project_436_2018.sav")
dt <- as.data.frame(dt)
dt <- na.omit(dt)
attach(dt)

#Answer to the question number 1.
t.test(dt$AGE[which(dt$GENDER=="M")], dt$AGE[which(dt$GENDER=="F")])
model3 <- aov(dt$AGE~dt$SCHOOL)
summary(model3)
#Answer to the question number 2.
install.packages("Hmisc")

library(Hmisc)
dt.cont <- data.frame(dt$AGE, dt$HEIGHT, dt$WEIGHT, dt$WAIST, dt$ARM)
pairs.default(dt.cont)
cor(dt.cont, method = "pearson")
cor_test <- rcorr(as.matrix(dt.cont, type= "pearson"))
cor_test$P
#Answer to the question number 3.

model2 <- lm(dt$bmi~dt$AGE+dt$HEIGHT+dt$WEIGHT+dt$WAIST+dt$ARM, data = dt)
install.packages("car")

library(car)

vif(model2)
library(MASS)

summary(model2)
step <- stepAIC(model2, direction = "backward")
step$anova

#Answer to the question number 3.
qqnorm(dt$bmi)
shapiro.test(dt$bmi)
model2 <- lm(dt$bmi~dt$AGE+dt$WAIST+dt$ARM)
summary(model2)
plot(model2)

# b<-boxcox(dt$bmi~1,lambda = seq(-2,0,1))
# b
# lamb<--1.25
# dt$bmib<-((dt$bmi^lamb)-1)/lamb
# 
# 
# model2 <- lm((1/(dt$bmi)^.5)~dt$AGE+dt$WAIST+dt$ARM)
# model2 <- lm(dt$bmib~dt$AGE+dt$WAIST+dt$ARM)
log(2.5)
plot(model2)


library(MASS)
bc <- boxcox(model2, lambda = seq(-0.25, 0.25))

p2 <- predict(model2, dt.set2, type = "response")

summary(model2)

model2$coefficients
coff <- as.vector(model2$coefficients)
n.coff <- (1/(coff))^0.5
n.coff




#Answer to the question number 4.
dt$bmi <- (dt$WEIGHT/(dt$HEIGHT/100)^2)
m <- median(dt$bmi)

dt$bmi_cat[dt$bmi < m] <- "Under Median"
dt$bmi_cat[dt$bmi >= m] <- "Over Median"
dt$over_median[dt$bmi_cat == "Under Median"] <- 0
dt$over_median[dt$bmi_cat == "Over Median"] <- 1
dt$over_median <- as.factor(dt$over_median)
model <- glm(formula = dt$over_median~ dt$AGE+ dt$WAIST +dt$ARM, family = "binomial")
summary(model)
dt.set <- data.frame(dt$AGE, dt$WAIST, dt$ARM)
prob <- predict(model, dt.set, type = "response")
dt$over_median_predict[prob < 0.5] <- 0
dt$over_median_predict[prob >= 0.5] <- 1
dt$over_median_predict <- as.factor(dt$over_median_predict)
table(dt$over_median, dt$over_median_predict)
#______________Yousuf End_______________

clin<-read.csv(file= "E:/222/CSE490/clinvar-conflicting/clinvar_conflicting.csv", header = T)
head(clin)
View(clin)
cl<-clin
create_report(cl)
colnames(cl)
View(cl)
cl<-cl[c(1:63246),]

txt2csv <- function(mydir, mycsvfilename){
  
  starting_dir <- getwd()
  
  # Get the names of all the txt files (and only txt files)
  myfiles <- list.files(mydir, full.names = TRUE, pattern = "*.txt")
  
  # Read the actual contexts of the text files into R and rearrange a little.
  
  # create a list of dataframes containing the text
  mytxts <- lapply(myfiles, readLines)
  
  # combine the rows of each dataframe to make one
  # long character vector where each item in the vector
  # is a single text file
  mytxts1lines <- unlist(mytxts)
  
  # make a dataframe with the file names and texts
  mytxtsdf <- data.frame(filename = basename(myfiles), # just use filename as text identifier
                         fulltext = mytxts1lines) # full text character vectors in col 2
  
  # Now write them all into a single CSV file, one txt file per row
  
  setwd(mydir) # make sure the CSV goes into the dir where the txt files are
  # write the CSV file...
  write.table(mytxtsdf, file = paste0(mycsvfilename, ".csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  # now check your folder to see the csv file
  message(paste0("Your CSV file is called ", paste0(mycsvfilename, ".csv"), ' and can be found in ', getwd()))
  
  # return original working directory
  setwd(starting_dir)
}

ctx<-read.csv(file = "E:/rough/Book1.csv",header = T)
View(ctx)
fix(ctx)
dir.create("E:/rough/naimul/ht/")
paste("Naimul", "Hasan", sep = " ,")

as.String(ctx[1,1])

gsub("[0-9]{1,4} (.*)","\\1",  as.String(ctx[1,3]))
ctxx<-as.String(c("Naimul","Hasan"))
removePunctuation(removeNumbers(paste(as.String(ctx[1,1]), collapse = " ")))

removePunctuation(removeNumbers(as.String(ctx[1,3])))
library(stringr)
c[1]<-removePunctuation(removeNumbers(paste(as.String(ctx[1,1]), collapse = " ")))
c[3]<-removePunctuation(removeNumbers(paste(as.String(ctx[1,3]), collapse = " ")))
c[2]<-removePunctuation(removeNumbers(paste(as.String(ctx[1,2]), collapse = " ")))

str_c(removePunctuation(removeNumbers(paste(as.String(ctx[1,1]), collapse = " "))),removePunctuation(removeNumbers(paste(as.String(ctx[1,3]), collapse = " ")))
, collapse = " ")


tm_map(ctx[1,],removeNumbers)
for(i in 1:2){
  for(j in 1:3)
  {
    c[j]<-removePunctuation(removeNumbers(paste(as.String(ctx[i,j]), collapse = " ")))
    
  }
  writeLines(c, 
             paste0("E:/rough/naimul/ht/outfile-", i, ".txt"))
}


txt <- c("here is", "some text", "to test", "this function with", "'including a leading quote", '"and another leading quote')
# make text files
dir.create("E:/rough/naimul/testdir")
getwd("E:/rough/naimul/testdir")
for(i in 1:length(txt)){
  writeLines(txt[i], paste0("E:/rough/naimul/testdir/outfile-", i, ".txt"))
}

# run the function and then look in the CSV file that is produced.
txt2csv("E:/rough/naimul/testdir", "theoutfile")

library(igraph)
install.packages("sna")

library(sna)

lin<-read.csv("E:/link.csv", header = T)
fix(lin)
li<-lin
rm(lin)

table(li$Student)
li$CourseName[li$Student==1]
tabulate(li[,c(4:10)])
s<-c(1,2,3)

for(i in 4:52)
{
  s[i]<-sum(li[,i])

}
s[7]
li[,7]
colnames(li[,7])
c<-colnames(li)
cs<-data.frame(c,s)
sort(cs[,2],decreasing = F)
View(cs)
l<-c(7,8,10,11,5,14)

hist(li$Viewers)
sort(li$Viewers, decreasing = T)

print(c(li$CourseName, sort(li$Viewers, decreasing = T)))


summary(s)

cat(sprintf("<set Value=\"%f\" Name=\"%s\" ></set>\n",sort(li$Viewers, decreasing = T),li$CourseName))

sli<-li[order(li$Viewers,decreasing = T),]
sli[c(1:10),c(2,3)]

#Warnings
assign("last.warning", NULL, envir = baseenv())
for(i in 1:107)
{
  if(1<=i&&i<=5)
  {
    dt1[,i]<-dt1[,i]
  }
  if(6<=i&&i<=107)
  {
    dt1[,i]<-as.numeric(unlist(dt1[,i]))
  }
}
#Download_Package
pkg.list = available.packages()
download.packages(pkgs = pkg.list, destdir = "F:/R_Packages")
rm(pkg.list)
install.packages("C:/Users/Staff Asia/Desktop/Rcpp_0.12.19.zip", repos=NULL)




#Small_Area_Estimation
sae<-read.csv("E:/222/Project Small Area Estimation/work.csv",header = T)
View(sae)
sae<-sae[c(1:58),]

length(unique(sae[,1]))
m.sam<-length(unique(sae[,1]))
y.sam<-sae[,2]
nd<-sae[,3]
sigma_square.i<-sae[,4]
x.sam<-sae[,5:6]
Xs.sam<-sae[,c(1,5,6)]
k<-c(dim(Xs.sam))[2]
fit1<-lm(y.sam~x.sam$HH_F+x.sam$HH_SIZE)
eibar<-fit1$residuals
l<-lm.influence(fit1)
l$hat
hi<-1-l$hat
comp2<-sum(sigma_square.i*hi)
delta<-comp1-comp2

matrix(0,3,3)


#___________Sentiment_Analysis_______________
#___________1&2________________

text <- c("Because I could not stop for Death -",
          "He kindly stopped for me -",
          "The Carriage held but just Ourselves -",
          "and Immortality")

text

library(dplyr)
text_df <- data_frame(line = 1:4, text = text)

text_df

text_df %>%
  unnest_tokens(word, text)

library(janeaustenr)
library(dplyr)
library(stringr)

original_books <- austen_books() %>%
  group_by(book) %>%
  mutate(linenumber = row_number(),
         chapter = cumsum(str_detect(text, regex("^chapter [\\divxlc]",
                                                 ignore_case = TRUE)))) %>%
  ungroup()

original_books

tidy_books <- original_books %>%
  unnest_tokens(word, text)

tidy_books

data(stop_words)

tidy_books <- tidy_books %>%
  anti_join(stop_words)
data(stop_words)


tidy_books %>%
  count(word, sort = TRUE) 

library(ggplot2)

tidy_books %>%
  count(word, sort = TRUE) %>%
  filter(n > 600) %>%
  mutate(word = reorder(word, n)) %>%
  ggplot(aes(word, n)) +
  geom_col() +
  xlab(NULL) +
  coord_flip()

install.packages("gutenbergr")

library(gutenbergr)

hgwells <- gutenberg_download(c(35, 36, 5230, 159))

tidy_hgwells <- hgwells %>%
  unnest_tokens(word, text) %>%
  anti_join(stop_words)

tidy_hgwells %>%
  count(word, sort = TRUE)

bronte <- gutenberg_download(c(1260, 768, 969, 9182, 767))

tidy_bronte <- bronte %>%
  unnest_tokens(word, text) %>%
  anti_join(stop_words)

tidy_bronte %>%
  count(word, sort = TRUE)

library(tidyr)
frequency <- bind_rows(mutate(tidy_bronte, author = "Bront� Sisters"),
                       mutate(tidy_hgwells, author = "H.G. Wells"), 
                       mutate(tidy_books, author = "Jane Austen")) %>% 
  mutate(word = str_extract(word, "[a-z']+")) %>%
  count(author, word) %>%
  group_by(author) %>%
  mutate(proportion = n / sum(n)) %>% 
  select(-n) %>% 
  spread(author, proportion) %>% 
  gather(author, proportion, `Bront� Sisters`:`H.G. Wells`)

library(scales)

# expect a warning about rows with missing values being removed
ggplot(frequency, aes(x = proportion, y = `Jane Austen`, color = abs(`Jane Austen` - proportion))) +
  geom_abline(color = "gray40", lty = 2) +
  geom_jitter(alpha = 0.1, size = 2.5, width = 0.3, height = 0.3) +
  geom_text(aes(label = word), check_overlap = TRUE, vjust = 1.5) +
  scale_x_log10(labels = percent_format()) +
  scale_y_log10(labels = percent_format()) +
  scale_color_gradient(limits = c(0, 0.001), low = "darkslategray4", high = "gray75") +
  facet_wrap(~author, ncol = 2) +
  theme(legend.position="none") +
  labs(y = "Jane Austen", x = NULL)


cor.test(data = frequency[frequency$author == "Bront� Sisters",],
         ~ proportion + `Jane Austen`)

cor.test(data = frequency[frequency$author == "H.G. Wells",], 
         ~ proportion + `Jane Austen`)


install.packages("tidytext")
library(tidytext)

sentiments

get_sentiments("afinn")
get_sentiments("bing")
get_sentiments("nrc")


library(janeaustenr)
library(dplyr)
library(stringr)

tidy_books <- austen_books() %>%
  group_by(book) %>%
  mutate(linenumber = row_number(),
         chapter = cumsum(str_detect(text, regex("^chapter [\\divxlc]", 
                                                 ignore_case = TRUE)))) %>%
  ungroup() %>%
  unnest_tokens(word, text)


nrc_joy <- get_sentiments("nrc") %>% 
  filter(sentiment == "joy")

tidy_books %>%
  filter(book == "Emma") %>%
  inner_join(nrc_joy) %>%
  count(word, sort = TRUE)

library(tidyr)

jane_austen_sentiment <- tidy_books %>%
  inner_join(get_sentiments("bing")) %>%
  count(book, index = linenumber %/% 80, sentiment) %>%
  spread(sentiment, n, fill = 0) %>%
  mutate(sentiment = positive - negative)

library(ggplot2)

ggplot(jane_austen_sentiment, aes(index, sentiment, fill = book)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~book, ncol = 2, scales = "free_x")

pride_prejudice <- tidy_books %>% 
  filter(book == "Pride & Prejudice")

pride_prejudice

afinn <- pride_prejudice %>% 
  inner_join(get_sentiments("afinn")) %>% 
  group_by(index = linenumber %/% 80) %>% 
  summarise(sentiment = sum(score)) %>% 
  mutate(method = "AFINN")

bing_and_nrc <- bind_rows(pride_prejudice %>% 
                            inner_join(get_sentiments("bing")) %>%
                            mutate(method = "Bing et al."),
                          pride_prejudice %>% 
                            inner_join(get_sentiments("nrc") %>% 
                                         filter(sentiment %in% c("positive", 
                                                                 "negative"))) %>%
                            mutate(method = "NRC")) %>%
  count(method, index = linenumber %/% 80, sentiment) %>%
  spread(sentiment, n, fill = 0) %>%
  mutate(sentiment = positive - negative)

bind_rows(afinn, 
          bing_and_nrc) %>%
  ggplot(aes(index, sentiment, fill = method)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~method, ncol = 1, scales = "free_y")


get_sentiments("nrc") %>% 
  filter(sentiment %in% c("positive", 
                          "negative")) %>% 
  count(sentiment)

get_sentiments("bing") %>% 
  count(sentiment)

bing_word_counts <- tidy_books %>%
  inner_join(get_sentiments("bing")) %>%
  count(word, sentiment, sort = TRUE) %>%
  ungroup()

bing_word_counts

library(ggplot2)

bing_word_counts %>%
  group_by(sentiment) %>%
  top_n(10) %>%
  ungroup() %>%
  mutate(word = reorder(word, n)) %>%
  ggplot(aes(word, n, fill = sentiment)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~sentiment, scales = "free_y") +
  labs(y = "Contribution to sentiment",
       x = NULL) +
  coord_flip()

custom_stop_words <- bind_rows(data_frame(word = c("miss"), 
                                          lexicon = c("custom")), 
                               stop_words)

custom_stop_words


library(wordcloud)

tidy_books %>%
  anti_join(stop_words) %>%
  count(word) %>%
  with(wordcloud(word, n, max.words = 100))

library(reshape2)

tidy_books %>%
  inner_join(get_sentiments("bing")) %>%
  count(word, sentiment, sort = TRUE) %>%
  acast(word ~ sentiment, value.var = "n", fill = 0) %>%
  comparison.cloud(colors = c("gray20", "gray80"),
                   max.words = 100)


PandP_sentences <- data_frame(text = prideprejudice) %>% 
  unnest_tokens(sentence, text, token = "sentences")

PandP_sentences$sentence[2]

austen_chapters <- austen_books() %>%
  group_by(book) %>%
  unnest_tokens(chapter, text, token = "regex", 
                pattern = "Chapter|CHAPTER [\\dIVXLC]") %>%
  ungroup()

austen_chapters %>% 
  group_by(book) %>% 
  summarise(chapters = n())

bingnegative <- get_sentiments("bing") %>% 
  filter(sentiment == "negative")

wordcounts <- tidy_books %>%
  group_by(book, chapter) %>%
  summarize(words = n())

tidy_books %>%
  semi_join(bingnegative) %>%
  group_by(book, chapter) %>%
  summarize(negativewords = n()) %>%
  left_join(wordcounts, by = c("book", "chapter")) %>%
  mutate(ratio = negativewords/words) %>%
  filter(chapter != 0) %>%
  top_n(1) %>%
  ungroup()

#_______________3_________________

library(dplyr)
library(janeaustenr)
library(tidytext)

book_words <- austen_books() %>%
  unnest_tokens(word, text) %>%
  count(book, word, sort = TRUE) %>%
  ungroup()

total_words <- book_words %>% 
  group_by(book) %>% 
  summarize(total = sum(n))

book_words <- left_join(book_words, total_words)

book_words

library(ggplot2)

ggplot(book_words, aes(n/total, fill = book)) +
  geom_histogram(show.legend = FALSE) +
  xlim(NA, 0.0009) +
  facet_wrap(~book, ncol = 2, scales = "free_y")

freq_by_rank <- book_words %>% 
  group_by(book) %>% 
  mutate(rank = row_number(), 
         `term frequency` = n/total)

freq_by_rank

freq_by_rank %>% 
  ggplot(aes(rank, `term frequency`, color = book)) + 
  geom_line(size = 1.1, alpha = 0.8, show.legend = FALSE) + 
  scale_x_log10() +
  scale_y_log10()


rank_subset <- freq_by_rank %>% 
  filter(rank < 500,
         rank > 10)

lm(log10(`term frequency`) ~ log10(rank), data = rank_subset)

freq_by_rank %>% 
  ggplot(aes(rank, `term frequency`, color = book)) + 
  geom_abline(intercept = -0.62, slope = -1.1, color = "gray50", linetype = 2) +
  geom_line(size = 1.1, alpha = 0.8, show.legend = FALSE) + 
  scale_x_log10() +
  scale_y_log10()

book_words <- book_words %>%
  bind_tf_idf(word, book, n)
book_words

book_words %>%
  select(-total) %>%
  arrange(desc(tf_idf))

book_words %>%
  arrange(desc(tf_idf)) %>%
  mutate(word = factor(word, levels = rev(unique(word)))) %>% 
  group_by(book) %>% 
  top_n(15) %>% 
  ungroup %>%
  ggplot(aes(word, tf_idf, fill = book)) +
  geom_col(show.legend = FALSE) +
  labs(x = NULL, y = "tf-idf") +
  facet_wrap(~book, ncol = 2, scales = "free") +
  coord_flip()

library(gutenbergr)
physics <- gutenberg_download(c(37729, 14725, 13476, 5001), 
                              meta_fields = "author")

install.packages("zoo")

library(zoo)

library(openxlsx)
emailID<-read.xlsx("E:/222/Our_Sites/One_ID.xlsx")

View(emailID)
emailID$user_registered
as.Date(emailID$user_registered, origin = "1899-12-30")
emailID$user_registered<-as.Date(emailID$user_registered, origin = "1899-12-30")
View(emailID)

IDandData<-read.xlsx("E:/222/Our_Sites/One_Courses.xlsx")
View(IDandData)
IDandData$date_recorded<-as.Date(IDandData$date_recorded, origin = "1899-12-30")
df <- merge(emailID, IDandData, by= c("user_id"), all.x= TRUE, all.y= TRUE)
View(df)

write.xlsx(df,"E:/222/Our_Sites/Combined_One.xlsx")

