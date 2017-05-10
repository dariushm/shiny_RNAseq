library(shiny)
library(monocle)

server = function(input, output, session) {
  options(shiny.maxRequestSize = 100 * 1024 ^ 2)
  
  #### Preprocessing I: Data upload
  observe({
    ma <- reactive({
      my_matrix <- input$matrix
      if (is.null(my_matrix)) {
        return(NULL)
      }
      read.csv(my_matrix$datapath, header = TRUE)
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
    
    ?reactiveValues()
    
    renamed_ph = reactive({
      pheno = ph()
      cols = 1:8
      cn_p = paste0('pheno', cols)
      print(names(input))
      cnames_p = lapply(cn_p, function(x) {
        input[[x]]
      })
      print(cnames_p)
      cnames_p = sapply(cnames_p, function(x) {
        if (is.null(x)) {
          return("")
        }
        x
      })
      print(cnames_p)
      if (all(cnames_p != "")) {
        colnames(pheno) = cnames_p
      }
      pheno
    })
    output$pd <- renderTable({
      pheno = renamed_ph()
      # print(pheno)
      if (is.null(pheno)) {
        p = NULL
      } else {
        p = head(pheno)
      }
      p
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
    observeEvent(input$acceptPheno, {
      updateTabsetPanel(session, "tabs", selected = "featuretab")
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
    
    renamed_fd = reactive({
      feature = fd()
      cols = 1:5
      cn_f = paste0('feature', cols)
      print(names(input))
      cnames_f = lapply(cn_f, function(x) {
        input[[x]]
      })
      print(cnames_f)
      cnames = sapply(cnames_f, function(x) {
        if (is.null(x)) {
          return("")
        }
        x
      })
      print(cnames_f)
      if (all(cnames_f != "")) {
        colnames(feature) = cnames_f
      }
      feature
    })
    output$fd <- renderTable({
      feature = renamed_fd()
      # print(pheno)
      if (is.null(feature)) {
        f = NULL
      } else {
        f = head(feature)
      }
      f
    })
    
    
    lapply(1:5, function(i) {
      output[[paste0('feature_col', i)]] <- renderUI({
        a = c(
          "X",
          "GeneShortName",
          "BioType",
          "#Cell",
          "UseForOrder"
        )
        df = fd()
        items <- colnames(df)
        selectInput(paste0('feature', i), a[i], c(items[i], items[-i]))
      })
    })
    
    
    observeEvent(input$acceptFeature, {
      updateTabsetPanel(session, "tabs", selected = "finalchecktab")
    })
    #cds <- newCellDataSet(as.matrix(ma),
    #                   phenoData = pheno,
    #                   featureData = feature)
    
    #output$selected=renderTable(head(cds))
    
    
  })
}
