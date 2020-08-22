setwd("C:/Users/Dolphin-Pc/Documents/R/win-library/3.4/knitr/shiny")
install.packages("shiny")
library(shiny)
shinyServer(
  pageWithSidebar(
    headerPanel("My First Shiny App"),
    
    sidebarPanel("Side Bar"),
    
    mainPanel("Main Panel")
  )
)

