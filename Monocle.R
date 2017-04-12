library(shiny)
ui = fluidPage(
  titlePanel("Monocle"),
  sidebarLayout(
    sidebarPanel(
      fileInput("matrix", "Expression Matrix", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      tags$hr(),
      checkboxInput('header', 'Header', TRUE),
      fileInput("phenodata", "Pheno Data", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      tags$hr(),
      checkboxInput('header', 'Header', TRUE),
      fileInput("featuredata", "Feature Data", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      tags$hr(),
      checkboxInput('header', 'Header', TRUE)
    ),
    
    mainPanel(
      tableOutput("pd"),
      tableOutput("fd")
    )
  )
)

server = function(input, output){
  options(shiny.maxRequestSize=100*1024^2) 
  ma <- reactive({
    matrix <- input$matrix
    if (is.null(matrix)) {return(NULL)}
    read.csv(matrix$datapath, header = TRUE)
  })
  
  ph <- reactive({
    pheno <- input$phenodata
    if (is.null(pheno)) {return(NULL)}
    read.csv(pheno$datapath, header = TRUE)
    
  })
  
  fe <- reactive({
    feature <- input$featuredata
    if (is.null(feature)) {return(NULL)}
    read.csv(feature$datapath, header = TRUE)
  })
  
  output$pd <- renderTable({
    pheno = ph()
    p = head(pheno)
  })
  
  output$fd <- renderTable({
    feature = fe()
    p = head(feature)
  })
  
}

shinyApp(ui = ui, server = server)
