library(shiny)

shinyUI({
  navbarPage("Monocle Dashboard",
             navbarMenu("Preprocess",
                        # Use example data
                        tabPanel("Upload data",
                                 sidebarLayout(
                                   sidebarPanel(
                                     # input matrix file
                                     fileInput("matrix", "Expression Matrix", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                     tags$hr(),
                                     checkboxInput('header', 'Header', TRUE),

                                     # input phenodata and modify the titles
                                     fileInput("phenodata", "Pheno Data", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                     tags$hr(),
                                     checkboxInput('header', 'Header', TRUE),

                                     # input featuredata and modify the titles
                                     fileInput("featuredata", "Feature Data", accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                     tags$hr(),
                                     checkboxInput('header', 'Header', TRUE)
                                     ),
                                   mainPanel(
                                     tabsetPanel(
                                       tabPanel("Pheno",
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
                                       ),
                                       tabPanel("Feature",
                                                splitLayout(
                                                  conditionalPanel(condition = "output.fdload", uiOutput(paste0('feature_col', 1)), div(style = "height:200px;")),
                                                  conditionalPanel(condition = "output.fdload", uiOutput(paste0('feature_col', 2)), div(style = "height:200px;")),
                                                  conditionalPanel(condition = "output.fdload", uiOutput(paste0('feature_col', 3)), div(style = "height:200px;")),
                                                  conditionalPanel(condition = "output.fdload", uiOutput(paste0('feature_col', 4)), div(style = "height:200px;")),
                                                  conditionalPanel(condition = "output.fdload", uiOutput(paste0('feature_col', 5)), div(style = "height:200px;")),
                                                  conditionalPanel(condition = "output.fdload", uiOutput(paste0('feature_col', 6)), div(style = "height:200px;")),
                                                  conditionalPanel(condition = "output.fdload", uiOutput(paste0('feature_col', 7)), div(style = "height:200px;")),
                                                  conditionalPanel(condition = "output.fdload", uiOutput(paste0('feature_col', 8)), div(style = "height:200px;"))
                                                ),
                                                tableOutput("fd")
                                       )
                                     )
                                   )
                                   )
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
