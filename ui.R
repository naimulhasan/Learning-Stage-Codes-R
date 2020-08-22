library(shiny)
library(shinydashboard)
shinyUI(
  dashboardPage(
    dashboardHeader(title="Histogram Site"),
    dashboardSidebar(
      sliderInput("bins","Number of Breaks",1,100,50), 
      menuItem("Dashboard"),
        menuSubItem("Dashboard Finance"),
        menuSubItem("Dashboard Sales"),
      menuItem("Detailed Analysis"),
      menuItem("Raw Data")
      ),
    dashboardBody(
      fluidRow(
        box(plotOutput("histogram"))
        #box(sliderInput("bins","Number of Breaks",1,100,50))  
        #checkboxInput()
  #numericInput()
  #textInput()
      )
    )
  )

#  ),
# runExample(),
#runExample("03_reactivity")
)

