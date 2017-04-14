library(shiny)
library(shinyBS)

shinyUI({
  navbarPage("Monocle Dashboard",
             navbarMenu("Preprocess",
                        # Use example data
                        tabPanel("Upload data",
                                 wellPanel("Use example data",
                                           selectInput("cdsExample", "Available datasets", 
                                                       c("hsmm", "lung"))),
                                 wellPanel(tableOutput("selected"))
                        ) 
             )
  )
})