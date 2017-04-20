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
             ),
             navbarMenu("Filtering", 
                        tabPanel("Select variable",
                                 wellPanel("XX",
                                           checkboxGroupInput("layer", "Select level to work with", 
                                                       c("pData", "fData", "eLevel")),
                                           selectInput("labels", "Rowheaders", 
                                                       c("pData", "fData", "eLevel"))))
                                 
  )
  )
})