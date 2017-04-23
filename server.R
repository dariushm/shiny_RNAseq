library(shiny)
<<<<<<< HEAD
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
=======
function(input, output, session) {
  options(shiny.maxRequestSize = 100*1024^2) 
  library(monocle)
  library(HSMMSingleCell)
  
  #### Preprocessing I: Data upload
  input$cdsExample=reactive({
    if (input$cdsExample == "hsmm") {
      data("HSMM_expr_matrix", "HSMM_gene_annotation", "HSMM_sample_sheet")
      cds <- newCellDataSet(as.matrix(HSMM_expr_matrix), 
                            new("AnnotatedDataFrame", data = HSMM_sample_sheet),
                            new("AnnotatedDataFrame", data = HSMM_gene_annotation),
                            expressionFamily = negbinomial())
      cds <- estimateSizeFactors(cds)
      cds <- estimateDispersions(cds)
    }
    if (input$cdsExample=="lung") {
          extPath=file.path(system.file(package = "monocle"), "extdata")
          load(file.path(extPath, "lung_phenotype_data.RData"))
           load(file.path(extPath, "lung_exprs_data.RData"))
           load(file.path(extPath, "lung_feature_data.RData"))
           cds <- newCellDataSet(lung_exprs_data[ , rownames(lung_phenotype_data)], 
                              new("AnnotatedDataFrame", data = lung_phenotype_data),
                              new("AnnotatedDataFrame", data = lung_feature_data),
                              lowerDetectionLimit=1,
                              expressionFamily = negbinomial.size())
           cds <- estimateSizeFactors(cds)
           cds <- estimateDispersions(cds)
    }
  })
  pRaw <- reactive({
    pheno <- input$phenoUser
    if (is.null(pheno)) {return(NULL)}
    read.csv(pheno$datapath, header = TRUE)
    })
  fRaw <- reactive({
    feature <- input$featureUser
    if (is.null(feature)) {return(NULL)}
    read.csv(feature$datapath, header = TRUE)
    })
  eRaw <- reactive({
    exprs <- input$exprUser
    if (is.null(matrix)) {return(NULL)}
    read.csv(matrix$datapath, header = TRUE)
    })
  
  ## Print top rows of loaded user files
  output$pd <- renderTable({
    pheno = ph()
    p = head(pheno)
    })
  output$fd <- renderTable({
    feature = fe()
    f = head(feature)
    })
  observeEvent(input$go.cds) {
    ## Put things in order
    eRaw <- eRaw[order(match(rownames(eRaw), rownames(fRaw))), 
                order(match(colnames(eRaw), rownames(pRaw)))]
  
    ## Check naming
    if (rownames(eRaw) != rownames(fRaw)) {
      badFeaturesWarning <- "Error: Gene names in the feature table do not match the gene names in the expression matrix"
      js_string <- 'alert("SOMETHING");'
      badFeaturesWarning <- sub("SOMETHING", badFeaturesWarning, js_string)
      session$sendCustomMessage(type = 'jsCode', list(value = badFeaturesWarning))
      }
    if (colnames(eRaw) != rownames(pRaw)) {
      badSamplesWarning <- "Error: Sample names in phenotype table do not match the sample names in the expression matrix"
      js_string <- 'alert("SOMETHING");'
      js_string <- sub("SOMETHING", badSamplesWarning, js_string)
      session$sendCustomMessage(type = 'jsCode', list(value = js_string))
      }
    pd <- new("AnnotatedDataFrame", data = pRaw)
    fd <- new("AnnotatedDataFrame", data = fRaw)
      
    ## Create CellDataSet object
    if (input$exprType == "countsRaw") {
      rawCtsWarning <- "Warning: Although Monocle can be used with raw read counts, these are not directly proportional to expression values unless you normalize them by length, so some Monocle functions could produce nonsense results. If you don't have UMI counts, We recommend you load up FPKM or TPM values instead of raw read counts."
      js_string <- 'alert("SOMETHING");'
      js_string <- sub("SOMETHING", rawCtsWarning, js_string)
      session$sendCustomMessage(type = 'jsCode', list(value = js_string))
      }
    if (input$exprType == "fpkm") {
      cds <- newCellDataSet(as.matrix(exprRaw), pd, fd, expressionFamily = tobit())
      }
    if (input$exprType == "fpkmLog") {
     cds <- newCellDataSet(as.matrix(exprRaw), pd, fd, expressionFamily = gaussianff())
      }
    if (input$exprType == "rpcSparse") {
      cds <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"), pd, fd, expressionFamily = negbinomial())
      }
    if (input$exprType == "rpc") {
      cds <- newCellDataSet(as.matrix(rpc_matrix), pd, fd, expressionFamily = negbinomial())
      }
    if (input$fpkmToRPC == TRUE && exprType == "fpkm" | "fpkmLog") {
      rpc_matrix <- relative2abs(cds)
      cds <- newCellDataSet(as(as.matrix(rpcData), "sparseMatrix"), pd, fd, expressionFamily = negbinomial())
      }
    if (input$fpkmToRPC == FALSE && exprType == "fpkm" | "fpkmLog") {
      fpkmWarning <- "Note: If you performed your single-cell RNA-Seq experiment using spike-in standards, you can convert these measurements into mRNAs per cell (RPC). RPC values are often easier to analyze than FPKM or TPM values, because have better statistical tools to model them. In fact, it's possible to convert FPKM or TPM values to RPC values even if there were no spike-in standards included in the experiment"
      js_string <- 'alert("SOMETHING");'
      js_string <- sub("SOMETHING", fpkmWarning, js_string)
      session$sendCustomMessage(type = 'jsCode', list(value = js_string))
    }
    
    ## Estimations
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    
    ## Setup cell attributes as filterable
    output$pdCols <- colnames(pData(cds))
    output$pdClass <- sapply(pData(cds), class)
  }
  
  #### Preprocessing II: Filtering/QC (iterative)
  ## Filter 
  observeEvent(input$go.distrib) {
    ## "Reset" -- need to have the removing genes/samples undo-able, but then set for going forward after user has ideal criteria
    fData(cds)$num_cells_expr <- NULL
    pData(cds)$num_genes_expr <- NULL
    pData(cds)$user_criteria <- NULL
    pData(cds)$Singlet <- NULL
    
    ## Filter by minimum expression
    fData(cds)$num_cells_expr <- Matrix::rowSums(exprs(cds) > input$min_expr)
    pData(cds)$num_genes_expr <- Matrix::colSums(exprs(cds) > input$min_expr)
    
    ## Filter by cell criteria 
    pData(cds)$user_criteria <- ifelse(
      pData(cds)$Cells.in.Well == 1 & 
        pData(cds)$Control == FALSE & 
        pData(cds)$Clump == FALSE, 
      TRUE, FALSE)
    ## Need to fix these to user input from ui.R
    
    ## Calculate total mRNA
    pData(cds)$Total_mRNAs <- Matrix::colSums(exprs(cds))
    upper_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) + 2*sd(log10(pData(cds)$Total_mRNAs)))
    lower_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) - 2*sd(log10(pData(cds)$Total_mRNAs)))
    
    ## Identify doublets
    pData(cds)$Singlet <- ifelse(pData(cds)$Total_mRNAs < input$doubletValue, TRUE, FALSE)
    
    ## Filter
    cdsFilt <- cds[fData(cds)$num_cells_expr == TRUE,
                   pData(cds)$num_genes_expr == TRUE & 
                     pData(cds)$user_criteria == TRUE &
                     pData(cds)$Singlet == TRUE]
    
    ## Plot distribution of total mRNA for each cell
    plot.distrib.cells <- ggplot(pData(cdsFilt), aes(x = Total_mRNAs, group = Hours, color = Hours)) +
      geom_density() +
      geom_vline(xintercept = lower_bound) +
      geom_vline(xintercept = upper_bound)
    renderPlot(plot.distrib.cells) 
    
    ## Plot distribution of fpkm for each gene
    exprLog <- log(exprs(cds[expressed_genes,]))
    melted_dens_df <- melt(Matrix::t(scale(Matrix::t(exprLog))))
    plot.distrib.genes <- ggplot(melted_dens_df, aes(x = value)) + 
      geom_density() +
      stat_function(fun = dnorm, size = 0.5, color = "red") +
      labs(x = "Standardized log(FPKM)", y = "Density")
    renderPlot(plot.distrib.genes)
    
    ## Add brush feature to remove outliers?
  }
  
  #### Analysis I: Classifying and counting cells
  ##### 1. Cell-type hierarchy (user provides genes that mark a population of cells)
  gene1_id <- reactive({
    rownames(subset(fData(cdsFilt), gene_short_name == input$gene1))
  })
  gene2_id <- reactive({
    rownames(subset(fData(cdsFilt), gene_short_name == input$gene2))
  })
  ## Need to print warning if gene name is not in fData(cds), and remind user of nomenclature
  cell1_id <- reactive({
    input$cell1
  })
  cell2_id <- reactive({
    input$cell2
  })
  
  ## Add cell type information
  cth <- newCellTypeHierarchy()
  cth <- addCellType(cth, cell1_id, classify_func = function (x) {
    x[gene1_id, ] >= 1
    })
  cth <- addCellType(cth, cell2_id, classify_func = function (x) {
    x[gene1_id, ] < 1 & x[gene2_id, ] > 1
    })
  
  ## Classify cells based on that information
  cds <- classifyCells(cds, cth, input$min_expr)
  
  ## Summarize results in table
  renderTable(table(pData(cds)$CellType))
  
  ## Summarize results with pie chart
  plot.supervised <- ggplot(pData(cds), aes(x = factor(1), fill = factor(CellType))) +
    geom_bar(width = 1) + 
    coord_polar(theta = "y") +
    theme(axis.title = element_blank())
  renderPlot(plot.supervised)

  ##### 2. Unsupervised cell clustering
  disp_table <- dispersionTable(cds)
  unsup_clustering_genes <- subset(disp_table, 
                                   mean_expression >= input$mean_expr & 
                                     dispersion_empirical >= 1 * dispersion_fit)
  
  ## Plot
  ## Need to add user input for what covariates/attributes to include in the model (eg Media)
  ## Need to make the rendering iterative/interactive
  cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
  plot_ordering_genes(cds)
  cds <- clusterCells(cds, num_clusters = input$ncell)
  plot_cell_trajectory(cds, 1, 2, color = "CellType")
  plot_cell_trajectory(cds, 1, 2, color = "Media")
  cds <- clusterCells(cds, residualModelFormulaStr = "~Media", num_clusters = input$ncell)
  plot_cell_trajectory(cds, 1, 2, color = "Media")
  cds <- clusterCells(cds, residualModelFormulaStr="~Media + num_genes_expr", num_clusters = input$ncell)
  plot_cell_trajectory(cds, 1, 2, color = "num_genes_expr")
  plot_cell_trajectory(cds, 1, 2, color = "CellType")

  #### 3. Semi-supervised cell clustering with known gene markers
  marker_diff <- markerDiffTable(cds[expressed_genes, ],
                                 cth,
                                 residualModelFormulaStr = "~Media", cores = 1)
  candidate_clustering_genes <- rownames(subset(marker_diff, qval < input$qval))
  marker_spec <- calculateMarkerSpecificity(cds[candidate_clustering_genes, ], cth)
  
  ## Print top 3 gene markers
  output$topGenes <- head(selectTopMarkers(marker_spec, 3))
  
  ## Use the top 20 markers for clustering
  semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 20)$gene_id)
  
  ## Plot
  cds <- setOrderingFilter(cds, semisup_clustering_genes)
  renderPlot(plot_ordering_genes(cds))
  cds <- clusterCells(cds,
                       num_clusters = input$ncell,
                       clustering_genes = semisup_clustering_genes,
                       residualModelFormulaStr = "~Media + num_genes_expr")
  renderPlot(plot_cell_trajectory(cds, 1, 2, 
                                  color = "CellType"))

  #### 4. Imputing cell type
  cds <- clusterCells(cds,
                       num_clusters = input$ncell,
                       frequency_thresh = input$freq,
                       cell_type_hierarchy = cth,
                       clustering_genes = rownames(subset(marker_diff, qval < qval)),
                       residualModelFormulaStr = "~Media + num_genes_expr")
  
  ## Plot
  renderPlot(plot_cell_trajectory(cds, 1, 2, 
                                  color = "CellType", 
                                  markers = c(gene1_id, gene2_id)))
  plot.imputed <- ggplot(pData(cds), aes(x = factor(1), fill = factor(CellType))) +
    geom_bar(width = 1) + coord_polar(theta = "y") +
    theme(axis.title = element_blank())
  renderPlot(plot.imputed)
  cds1 <- cds[ , pData(cds)$CellType == cds1]
  cds1 <- estimateDispersions(cds1)

  ##### Analysis II: Constructing single-cell trajectories ("Pseudotime").
  ##### 1. Unsupervised ordering
  mean_expr <- reactive({
    input$mean_expr
  })
  qval <- reactive({
    input$qval
  })
  num_cells_expressed <- reactive({
    input$num_cells_expressed
  })
  genesUser <- reactive({
    as.character(input$genesUser)
    ## Check to see if genes exist in fData!
  })
  diff_test_res <- differentialGeneTest(cds1[expressed_genes, ],
                                        fullModelFormulaStr = "~Media")
  ordering_genes <- rownames (subset(diff_test_res, qval < qval))
  disp_table <- dispersionTable(cds1)
  ordering_genes <- subset(disp_table,
                           mean_expression >= input$mean_exp &
                             dispersion_empirical >= 2 * dispersion_fit)$gene_id
  cds1 <- setOrderingFilter(cds1, ordering_genes)
  
  ## Plot
  renderPlot(plot_ordering_genes(cds1))
  cds1 <- reduceDimension(cds1, max_components = 2)
  cds1 <- orderCells(cds1, reverse = FALSE)
  
  ## Plot
  renderPlot(plot_cell_trajectory(cds1, color_by = "Hours"))
  
  ## Plot subset of genes
  cds1_expr_genes <- rownames(subset(fData(cds1), num_cells_expressed >= num_cells_expressed))
  cds1Filt <- cds1[cds1_expr_genes,]
  my_genes <- rownames(subset(fData(cds1Filt),
                               gene_short_name %in% genesUser))
  cds_subset <- cds1Filt[my_genes, ]
  renderPlot(plot_genes_in_pseudotime(cds_subset, color_by = "Hours"))
  
  #### 2. Selecting genes based on PCA loading
  expr_Filt <- t(t(exprs(cds1Filt)/pData(cds1Filt)$Size_Factor))
  nz_genes <- which(expr_Filt != 0)
  expr_Filt[nz_genes] <- log(expr_Filt[nz_genes] + 1)
  
  ## Calculate the variance across genes without converting to a dense matrix:
  expr_means <- Matrix::rowMeans(expr_Filt)
  expr_vars <- Matrix::rowMeans((expr_Filt - expression_means)^2)
  
  ## Filter out genes that are constant across all cells:
  genes_to_keep <- expr_vars > 0
  expr_filt <- expr_Filt[genes_to_keep, ]
  expr_means <- expr_means[genes_to_keep]
  expr_vars <- expr_vars[genes_to_keep]
  
  ## Take top PCA loading genes
  irlba_pca_res <- irlba(t(expr_Filt),
                         nu = 0,
                         center = expr_means,
                         scale = sqrt(expr_vars),
                         right_only = TRUE)$v
  rownames(irlba_pca_res) <- rownames(expr_Filt)
  PC2_genes <- names(sort(abs(irlba_pca_res[ , 2]), decreasing = TRUE))[1:100]
  PC3_genes <- names(sort(abs(irlba_pca_res[ , 3]), decreasing = TRUE))[1:100]
  ordering_genes <- union(PC2_genes, PC3_genes)
  cds1Filt <- setOrderingFilter(cds1Filt, ordering_genes)
  cds1Filt <- reduceDimension(cds1Filt, max_components = 2)
  cds1Filt <- orderCells(cds1Filt, reverse = FALSE)
  
  ## Plot
  renderPlot(plot_cell_trajectory(cds1Filt, color_by = "Hours"))
  
  
  #### 3. Semi-supervised ordering with known marker genes
  cth <- newCellTypeHierarchy()
  cth <- addCellType(cth, cell1_id, classify_func = function (x) {
    x[gene1_id, ] >= 1
    })
  cth <- addCellType(cth, cell2_id, classify_func = function (x) {
    x[gene2_id, ] >= 1
    })
  marker_diff <- markerDiffTable(cds1[expressed_genes, ],
                                 cth,
                                 cores = 1)
  semisup_clustering_genes <- rownames(subset(marker_diff, qval < qval))
  
  ## Print number
  renderTable(length(semisup_clustering_genes))
  
  cds1 <- setOrderingFilter(cds1, semisup_clustering_genes)
  renderPlot(plot_ordering_genes(cds1))
  cds1 <- reduceDimension(cds1, max_components = 2)
  cds1 <- orderCells(cds1, reverse = TRUE)
  
  ## Plot
  renderPlot(plot_cell_trajectory(cds1, color_by = "Hours"))
  
  ## Plot subset of genes
  cds1Filt <- cds1[expressed_genes, ]
  my_genes <- rownames(subset(fData(cds1Filt), gene_short_name %in% genesUser))
  cds_subset <- cds1Filt[my_genes, ]
  renderPlot(plot_genes_in_pseudotime(cds_subset, color_by = "Hours"))
  
  #### Reconstructing branched trajectories
  ## Need to adjust for user data
  lung <- load_lung()
  
  color_by <- reactive({
      input$color_by
    })
  facet_by <- reactive({
    input$facet_by
  })
  root_state <- reactive({
    input$nroot
  })
  genesUser <- reactive({
    as.character(input$genesUser)
    ## Check to see if genes exist in fData!
  })
  
  ## Plot
  ## Need to make this reactive, add facet as reactive
  plot.traject <- plot_cell_trajectory(orderCells(lung, root_state = root_state),
                                       color_by = color_by)
  renderPlot(plot.traject)
  renderPlot(plot.traject + facet_wrap(~facet_by))

  
  ##### Analysis III: Differential expression analysis.
  ##### Basic
  marker_genes <- rownames(subset(fData(csd), gene_short_name %in% genesUser))
  diff_test_res <- differentialGeneTest(cds1Filt[marker_genes, ],
                                        fullModelFormulaStr = "~Media")
  sig_genes <- subset(diff_test_res, qval < qval)
  
  ## Print results
  renderTable(sig_genes[ , c("gene_short_name", "pval", "qval")])
  
  ## Plot
  cds_subset <- cds1Filt[rownames(subset(fData(cds1Filt), gene_short_name %in% genesUser)), ]
  renderPlot(plot_genes_jitter(cds_subset, grouping = "Media", ncol = 2))

  #### Finding genes that distinguish cell type or state
  to_be_tested <- rownames(subset(fData(cds), gene_short_name %in% genesUser))
  cds_subset <- cds[to_be_tested, ]
  diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~CellType")
  
  ## Print results
  renderTable(diff_test_res[ , c("gene_short_name", "pval", "qval")])
  
  ## Plot
  renderPlot(plot_genes_jitter(cds_subset, 
                               grouping = "CellType", 
                               color_by = "CellType",
                               nrow = 1, 
                               ncol = NULL, 
                               plot_trend = TRUE))
  full_model_fits <- fitModel(cds_subset, modelFormulaStr = "~CellType")
  reduced_model_fits <- fitModel(cds_subset, modelFormulaStr = "~1")
  diff_test_res <- compareModels(full_model_fits, reduced_model_fits)
  
  ## Print results
  renderTable(diff_test_res)

  #### Finding genes that change as a function of pseudotime
  to_be_tested <- rownames(subset(fData(cds), gene_short_name %in% genesUser))
  cds_subset <- cds1[to_be_tested, ]
  diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)")
  
  ## Print results
  renderTable(diff_test_res[ , c("gene_short_name", "pval", "qval")])
  
  ## Plot
  renderPlot(plot_genes_in_pseudotime(cds_subset, color_by = "Hours"))
  
  ## Clustering genes by pseudotemporal expression pattern
  diff_test_res <- differentialGeneTest(cds1[marker_genes, ],
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
  sig_gene_names <- rownames(subset(diff_test_res, qval < qval))
  
  ## Plot
  plot.pseudoCluster <- plot_pseudotime_heatmap(cds1[sig_gene_names, ], 
                                                num_clusters = 2, 
                                                cores = 1, 
                                                show_rownames = TRUE)
  renderPlot(plot.pseudoCluster)
  
  
  #### Multi-factorial differential gene expression analysis
  to_be_tested <- rownames(subset(fData(cds), gene_short_name %in% genesUser))
  cds_subset <- cds[to_be_tested, ]
  diff_test_res <- differentialGeneTest(cds_subset,
                                        fullModelFormulaStr = "~CellType + Hours",
                                        reducedModelFormulaStr = "~Hours")
  diff_test_res[ , c("gene_short_name", "pval", "qval")]
  renderPlot(plot_genes_jitter(cds_subset,
                               grouping = "Hours", 
                               color_by = "CellType", 
                               plot_trend = TRUE) +
               acet_wrap(~feature_label, scales = "free_y"))

  #### Analyzing branches in single-cell trajectories
  ## Need to adjust for user data
  lung <- load_lung()
  
  renderPlot(plot_cell_trajectory(lung, 
                                  color_by = "Time"))
  BEAM_res <- BEAM(lung, branch_point = 1, cores = 1)
  BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
  BEAM_res <- BEAM_res[ , c("gene_short_name", "pval", "qval")]
  
  ## Plot heatmap
  renderPlot(plot_genes_branched_heatmap(lung[rownames(subset(BEAM_res, qval < qval)), ],
                                         branch_point = 1,
                                         num_clusters = 4,
                                         cores = 1,
                                         use_gene_short_name = TRUE,
                                         show_rownames = TRUE))
  lung_genes <- rownames(subset(fData(lung), gene_short_name %in% genesUser))
  
  ## Plot pseudotime
  renderPlot(plot_genes_branched_pseudotime(lung[lung_genes, ],
                                            branch_point = 1,
                                            color_by = "Time",
                                            ncol = 1))
}
  
>>>>>>> a7e119605885a745b37078957ba7ec736eba7c71
