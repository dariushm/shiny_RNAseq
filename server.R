library(shiny)
library(shinyBS)
library(monocle)
library(HSMMSingleCell)

function(input, output, session) {
  options(shiny.maxRequestSize = 100*1024^2) 
  
  #### Preprocessing I: Data upload
  observeEvent(input$cdsExample, {
    if (input$cdsExample == "hsmm") {
      data("HSMM_expr_matrix", "HSMM_gene_annotation", "HSMM_sample_sheet")
      cds <- newCellDataSet(as.matrix(HSMM_expr_matrix), 
                            new("AnnotatedDataFrame", data = HSMM_sample_sheet),
                            new("AnnotatedDataFrame", data = HSMM_gene_annotation),
                            lowerDetectionLimit = 0.1,
                            expressionFamily = negbinomial())
      cds <- estimateSizeFactors(cds)
      cds <- estimateDispersions(cds)
      
      output$selected=renderTable(head(pData(cds)))
    }
    if (input$cdsExample == "lung") {
      extPath=file.path(system.file(package="monocle"), "extdata")
      load(file.path(extPath, "lung_phenotype_data.RData"))
      load(file.path(extPath, "lung_exprs_data.RData"))
      load(file.path(extPath, "lung_feature_data.RData"))
      cds=newCellDataSet(lung_exprs_data[, rownames(lung_phenotype_data)], 
                         new("AnnotatedDataFrame", data = lung_phenotype_data),
                         new("AnnotatedDataFrame", data = lung_feature_data),
                         lowerDetectionLimit = 1,
                         expressionFamily=negbinomial.size())
      cds <- estimateSizeFactors(cds)
      cds <- estimateDispersions(cds)
      
      output$selected=renderTable(head(pData(cds)))
    }
  })
  
  ###Step2.Select data that should be filtered
  observeEvent(input$layer, {
    x=input$layer
    if(is.null(x))
      x=character(0)
    if (x=="fData") {
      y=colnames(fData(cds))
      updateSelectInput(session, "labels",
                        label=paste("Select input label", length(y)),
                        choices=y, 
                        selected=tail(y, 1))
      }
    if (x=="pData") {
      y=colnames(pData(cds))
      updateSelectInput(session, "labels",
                        label=paste("Select input label", length(y)),
                        choices=y, 
                        selected=tail(y, 1))
    }
    if (x=="eLevel") {
      y=rownames(exprs(cds))
      updateSelectInput(session, "labels",
                        label=paste("Select input label", length(y)),
                        choices=y, 
                        selected=tail(y, 1))
    }
    
  ####Step3. Brushed points for data selection

  })
}
