library(shiny)

# Define UI for dataset viewer application
shinyUI(fluidPage(
  
  # Application title.
  titlePanel("Mouse to Human Prediction Interface"),

  sidebarLayout(
    sidebarPanel(
      selectInput("topic", "Narrow Topic?",choices=tests2,selected="ALL")
      
    ),
  
    mainPanel(
      h2("Data for Phenotypes"),
      tableOutput("summary"),
      plotOutput("histogram")
      
    )
  )
))
