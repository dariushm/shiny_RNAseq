library(shiny)
ui = fluidPage(
  titlePanel("Monocle"),
  sidebarLayout(
    sidebarPanel(
      fileInput("phenodata", "Pheno Data", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
    tags$hr(),
    checkboxInput('header', 'Header', TRUE)
    ),
    mainPanel(
      splitLayout(
        conditionalPanel(condition = "output.phload", uiOutput(paste0('pheno_col', 1)), div(style = "height:200px;")),
        conditionalPanel(condition = "output.phload", uiOutput(paste0('pheno_col', 2)), div(style = "height:200px;")),
        conditionalPanel(condition = "output.phload", uiOutput(paste0('pheno_col', 3)), div(style = "height:200px;")),
        conditionalPanel(condition = "output.phload", uiOutput(paste0('pheno_col', 4)), div(style = "height:200px;")),
        conditionalPanel(condition = "output.phload", uiOutput(paste0('pheno_col', 5)), div(style = "height:200px;")),
        conditionalPanel(condition = "output.phload", uiOutput(paste0('pheno_col', 6)), div(style = "height:200px;")),
        conditionalPanel(condition = "output.phload", uiOutput(paste0('pheno_col', 7)), div(style = "height:200px;")),
        conditionalPanel(condition = "output.phload", uiOutput(paste0('pheno_col', 8)), div(style = "height:200px;"))
      ),
      tableOutput("pd")
      # lapply(1:8, function(i){
      #   conditionalPanel(condition = "output.phload", uiOutput(paste0('pheno_col', i)))
      # })
    )
  )
)


server = function(input, output){
  options(shiny.maxRequestSize=100*1024^2) 
  
  ph <- reactive({
    infile = input$phenodata
    if(is.null(infile)) return(NULL)
    read.csv(infile$datapath, header = TRUE)
  })
  
  output$phload <- reactive({
    return(!is.null(ph()))
  })
  outputOptions(output, 'phload', suspendWhenHidden=FALSE)
  
  output$pd <- renderTable({
    pheno = ph()
    p = head(pheno)
  })
  
  
  lapply(1:8, function(i) {
    output[[paste0('pheno_col', i)]] <- renderUI({
    a = c("X", "Library", "Well", "Hours", "Media", "Mapped.Fragments", "Pseudotime", "State")
      df = ph()
      items <- colnames(df)
      selectInput(paste0('pheno', i), a[i], c(items[i], items[-i]))
    })
  })
  
  
}


shinyApp(ui = ui, server = server)
