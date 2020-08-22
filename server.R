library(shiny)
library(shinydashboard)

shinyServer(function(input,output)
{
  output$histogram<-renderPlot({
   #hist(10*(read.table(file="https://mirrors.netix.net/sourceforge/i/ir/irisdss/IRIS.csv", header = T, sep=",")$X5.1),breaks = input$bins)
    hist(iris$Sepal.Length,breaks = input$bins)
      })
})