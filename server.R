library(shiny)
library(monocle)

server = function(input, output) {
  options(shiny.maxRequestSize = 100 * 1024 ^ 2)
  ma <- reactive({
    matrix <- input$matrix
    if (is.null(matrix)) {
      return(NULL)
    }
    read.csv(matrix$datapath, header = TRUE)
  })

  ph <- reactive({
    pheno <- input$phenodata
    if (is.null(pheno)) {
      return(NULL)
    }
    read.csv(pheno$datapath, header = TRUE)
  })

  output$phload <- reactive({
    return(!is.null(ph()))
  })
  outputOptions(output, 'phload', suspendWhenHidden = FALSE)

  output$pd <- renderTable({
    pheno = ph()
    p = head(pheno)
  })

  lapply(1:8, function(i) {
    output[[paste0('pheno_col', i)]] <- renderUI({
      a = c(
        "X",
        "Library",
        "Well",
        "Hours",
        "Media",
        "Mapped.Fragments",
        "Pseudotime",
        "State"
      )
      df = ph()
      items <- colnames(df)
      selectInput(paste0('pheno', i), a[i], c(items[i], items[-i]))
    })
  })

  fd <- reactive({
    feature <- input$featuredata
    if (is.null(feature)) {
      return(NULL)
    }
    read.csv(feature$datapath, header = TRUE)

  })

  output$fdload <- reactive({
    return(!is.null(fd()))
  })
  outputOptions(output, 'fdload', suspendWhenHidden = FALSE)

  output$fd <- renderTable({
    feature = fd()
    f = head(feature)
  })

  lapply(1:8, function(i) {
    output[[paste0('feature_col', i)]] <- renderUI({
      a = c(
        "X",
        "Library",
        "Well",
        "Hours",
        "Media",
        "Mapped.Fragments",
        "Pseudotime",
        "State"
      )
      df = fd()
      items <- colnames(df)
      selectInput(paste0('feature', i), a[i], c(items[i], items[-i]))
    })
  })

}
