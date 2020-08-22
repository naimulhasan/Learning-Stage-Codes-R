library(ggplot2)
mt<-data.frame(mtcars)
mt
q111<-qplot(mt$cyl, 
      geom = "bar",
      fill= I("Pink"),
      colour=I("Blue"),
      xlab = "Cylinders",
      ylab= "Number of Vehicles",
      main= "Cylinders in mtcars"
      )
qplot(mt$hp, 
      geom = "histogram",
      binwidth=25,
      colour=I("black"),
      xlim = c(50,250),
      xlab = "HORSEPOWER",
      ylab = "Number of Cars",
      #alpha=T(0)
     main = "Histogram"
)
#Pie Diagram
barp<-ggplot(mt,aes(x=1,y=sort(mt$carb),fill=sort(mt$carb))) +
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
    labs(y="Carburetors")
)
count <- table(mtcars$cyl)
barplot(count)
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}
library(ggplot2)
qplot(mtcars$cyl, geom = "bar")
ggplot(data= mtcars, aes(mtcars$hp))  + 
  geom_histogram(breaks=seq(40, 350, by = 25),  
                 colour = I("black"),
                 aes(fill=..count..))
summary(mtcars)
mtcars
qplot(mtcars$cyl, geom = "bar", fill = I("black"), colour = I("red"))
qplot(mtcars$cyl, geom = "bar", fill = I("blue"), xlab = "Cylinders", ylab = "Number of Vehicles")
qplot(mtcars$cyl, geom = "bar", fill = I("blue"), xlab = "Cylinders", ylab = "Number of Vehicles", main = "Cylinders in mtcars")
qplot(mtcars$hp, geom="histogram")
qplot(mtcars$hp, geqplot(mtcars$hp, geom="histogram", binwidth = 25, colour = I("black")))
qplot(mtcars$hp, geom="histogram", binwidth = 25, colour = I("black"), xlab = "Horsepower", ylab= "Number of Cars", alpha = I(0))
qplot(mtcars$hp, geom="histogram", binwidth = 25, colour = I("black"), xlab = "Horsepower", ylab= "Number of Cars", alpha = I(0),
     main = "Histogram")
barp <- ggplot(mtcars, aes(x=1, y=sort(mtcars$carb), fill=sort(mtcars$carb))) +
        geom_bar(stat="identity")
print(barp)
#Scatter Diagram
qplot(mpg, wt, data=mtcars)
ggplot(mtcars, aes(x = mpg, y = wt)) + geom_point(shape=1)
ggplot(mtcars,aes(x=mpg,y=wt)) + geom_point(shape=19)
ggplot(mtcars,aes(x=mpg,y=wt,shape = cyl)) + geom_point() + scale_shape_identity()
mtcars$cylFactor<- factor(mtcars$cyl)
ggplot(mtcars,aes(x=mpg,y=wt,shape = cylFactor)) + geom_point()
ggplot(mtcars,aes(x=mpg,y=wt)) + geom_point(shape=19, colour="blue")
ggplot(mtcars,aes(x=mpg,y=wt,color = cyl))+ geom_point(shape=19)
ggplot(mtcars,aes(x=mpg,y=wt,color = cylFactor)) + geom_point(shape=19)
ggplot(mtcars,aes(x=mpg,y=wt,color = cylFactor)) + geom_point(shape=19) + labs(colour = "Cylinders")
ggplot(mtcars,aes(x=mpg,y=wt,color = cylFactor)) + geom_point(shape=19) + xlab("Miles per Gallon ") + ylab("Weight") +
  labs(colour = "Cylinders") + ggtitle("Scatterplot")
#Line Diargam
EuStockDF <- as.data.frame(EuStockMarkets)
ggplot(EuStockDF,aes(x=c(1:nrow(EuStockDF)), y = DAX)) + geom_line()
ggplot(EuStockDF,aes(x=c(1:nrow(EuStockDF)), y = DAX)) + geom_line(size=1.5) + labs(x = "Stocks")
ggplot(EuStockDF,aes(x=c(1:nrow(EuStockDF)), y = DAX)) + geom_line(size=1.5, colour="light blue") + labs(x = "Time", y = "Stocks")

dax_smi_plot <- ggplot() +
  geom_line(data = EuStockDF,aes(x=c(1:nrow(EuStockDF)), y = DAX), size = 1.5, colour="light blue") +
  geom_line(data = EuStockDF,aes(x=c(1:nrow(EuStockDF)), y = SMI), size = 1.5, colour = "red") +
  labs(x = "Time", y = "Stocks")
print(dax_smi_plot)

all_stocks <- ggplot() +
  geom_line(data = EuStockDF,aes(x=c(1:nrow(EuStockDF)), y = DAX), size=1, colour="light blue") +
  geom_line(data = EuStockDF,aes(x=c(1:nrow(EuStockDF)), y = SMI), size =1, colour = "red") + 
  geom_line(data = EuStockDF,aes(x=c(1:nrow(EuStockDF)), y = CAC), size =1, colour = "purple") + 
  geom_line(data = EuStockDF,aes(x=c(1:nrow(EuStockDF)), y = FTSE), size =1, colour = "green") +
  labs(x = "Time", y = "Stocks")
print(all_stocks)

legend_stocks <- all_stocks + xlab("Days") + ylab("Price") + ggtitle("Eu Stocks")
print(legend_stocks)

ggplot(mtcars,aes(x=mpg,y=wt))  + geom_point(shape=19) +
  geom_smooth(method="lm", se= FALSE, color = "red")

#se = TRUE  -> confidence interval appear (default = true)
ggplot(mtcars,aes(x=mpg,y=wt,color = cylFactor)) + geom_point(shape=19) + geom_smooth(method="lm", se= TRUE, color = "red")
#Regression Line
ggplot(mtcars,aes(x=mpg,y=wt,color = cylFactor)) + geom_point(shape=19) + 
  geom_smooth(method="lm", se= TRUE, color = "red") + xlab("Miles per Gallon ") + 
  ylab("Weight") +  labs(colour = "Cylinders") + ggtitle("Linear Regression")

#Gaussian Regression line(Curve)
ggplot(mtcars,aes(x=mpg,y=wt,color = cylFactor)) + geom_point(shape=19) + 
  geom_smooth(method="auto", se= TRUE, color = "red") + xlab("Miles per Gallon ") + 
  ylab("Weight") +  labs(colour = "Cylinders") + ggtitle("Gaussian Regression")
#Word Cloud

install.packages("tm")
install.packages("wordcloud")
library(tm)
library(wordcloud)
?Corpus

dir.create("G:/Programming/R code/wordcloud")
download.file("https://ibm.box.com/shared/static/cmid70rpa7xe4ocitcga1bve7r0kqnia.txt",
              destfile = "G:/Programming/R code/wordcloud/churchill_speeches.txt", quiet = TRUE)
dirPath<-"G:/Programming/R code/wordcloud"
dirPath<-"E:/rough"

speech<-Corpus(DirSource(dirPath))
View(speech)
?tm_map

inspect(speech)
speech<-tm_map(speech,content_transformer(tolower))
speech<-tm_map(speech,removeNumbers)
speech<-tm_map(speech,removeWords,stopwords("english"))
#_________2
dir.create("G:/Programming/R code/wordcloud/2")
download.file("https://ibm.box.com/shared/static/cmid70rpa7xe4ocitcga1bve7r0kqnia.txt",
              destfile = "G:/Programming/R code/wordcloud/churchill_speeches.txt", quiet = TRUE)
dirPath<-"G:/Programming/R code/wordcloud/2"
speech<-Corpus(DirSource(dirPath))
inspect(speech)
speech<-tm_map(speech,content_transformer(tolower))
speech<-tm_map(speech,removeNumbers)
speech<-tm_map(speech,removeWords,stopwords("english"))
#speech<-tm_map(speech,removeWords,c("any","words i want"))
speech<-tm_map(speech,removePunctuation)
speech<-tm_map(speech,stripWhitespace)
dtm<-TermDocumentMatrix(speech)
m<-as.matrix(dtm)
v<-sort(rowSums(m),decreasing = T)
d<-data.frame(word=names(v),freq=v)
head(d)
wordcloud(words = d$word,freq = d$freq)
wordcloud(words = d$word,freq = d$freq, min.freq = 1,max.words = 250)
wordcloud(words = d$word,freq = d$freq, 
          min.freq = 1,max.words = 250,
          colors = brewer.pal(8,"Dark2")
          )
wordcloud(words = d$word,freq = d$freq, 
          min.freq = 1,max.words = 250,
          colors = brewer.pal(8,"Dark2"),
          random.order = F
)
#Radar Charts
library(ggplot2)
install.packages("devtools")
library(devtools)
devtools::install_github("ricardo-bion/ggradar", 
                         dependencies=T)
install.packages("ggradar")
install.packages("dplyr")
install.packages("scales")
library(ggradar)
library(dplyr)
library(scales)
mtcars
library(scales)

mtcars %>%
  add_rownames( var = "group" ) %>%
  mutate_each(funs(rescale), -group) %>%
  tail(4) %>% select(1:10) -> mtcars_radar
mtcars %>%
  add_rownames(var="group") %>%
  mutate_at(funs(rescale),-group) %>%
  head(3) %>% select(1:10) -> mtcars_radar
options(warn = -1)
ggradar(mtcars_radar)
IRkernel::set_plot_options(width=950,height=600,units='px')
ggradar(mtcars_radar)

#Waffle Chart
library(ggplot2)
install.packages("IRkernel")
install.packages("waffle")
library(waffle)
library(IRkernel)
expenses <- c(`Health ($43,212)`=43212, `Education ($113,412)`=113412,
              `Transportation ($20,231)`=20231, `Entertainment ($28,145)`=28145)

waffle(expenses/1235, rows=5, size=0.3, 
       colors=c("#c7d4b6", "#a3aabd", "#a0d0de", "#97b5cf"), 
       title="Imaginary Household Expenses Each Year", 
       xlab="1 square = $934")

IRkernel::set_plot_options(width=950, height=600, units='px')
waffle(expenses/1235, rows=5, size=0.3, 
       colors=c("#c7d4b6", "#a3aabd", "#a0d0de", "#97b5cf"), 
       title="Imaginary Household Expenses Each Year", 
       xlab="1 square = $934")
#Box Plot
set.seed(123)
set_a<-rnorm(200,mean=1,sd=2)
set_b<-rnorm(200,mean=0,sd=2)
df<-data.frame(label=factor(rep(c("A","B"),each=200)),value =c(set_a,set_b))
df
library(ggplot2)
install.packages("plotly")
library(plotly)
ggplot(df,aes(x=label,y=value))+geom_boxplot()
ggplotly()
qplot(factor(mtcars$cyl),mtcars$mpg,geom="boxplot")

cars<-ggplot(mtcars,aes(factor(cyl),mpg))
cars+geom_boxplot()
#Maps

install.packages("leaflet")
library(leaflet)
map<-leaflet()
map
map<-leaflet() %>% addTiles()
map<-leaflet() %>% addTiles() %>% addMarkers(lng=-73.9851,lat=40.7589)
map<-leaflet() %>% addTiles() %>% 
  addMarkers(lng=-73.9851,lat=40.7589,popup='Times Square')
map<-leaflet() %>% addProviderTiles("Stamen.Watercolor") %>%
        addMarkers(lng=2.2945,lat=48.8584,
                   popup="Eiffel Tower")
#Stamen.TonerHybrid
#https://leaflet-extras.github.io/leaflet-providers/preview/
head(quakes)
map<-leaflet(quakes) %>% addTiles() %>%
      addCircleMarkers(lng=quakes$long,lat=quakes$lat)
map<-leaflet(quakes) %>% addTiles() %>%
  addMarkers(lng=quakes$long,lat=quakes$lat)
map<-leaflet(quakes) %>% addTiles() %>%
  addCircles(lng=quakes$long,lat=quakes$lat)
install.packages("htmlwidgets")
library(htmlwidgets)
install.packages("IRdisplay")
library(IRdisplay)

map<-leaflet(quakes) %>% addTiles() %>%
  addMarkers(clusterOptions=markerClusterOptions())
saveWidget(map,file="map7.html",selfcontained=F)
display_html(paste("<iframe src=' ",'map7.html',"'width='100%' height='300'","/>"))

map<-leaflet(quakes) %>% addTiles() %>%
  addMarkers(lng=86.92,lat=27.99,
  popup ="Mount Everest") %>%
  addRectangles(86.9,27.95,87,28.05)


png('Graph.png', 1000, 1000, type='cairo')
dev.off()







