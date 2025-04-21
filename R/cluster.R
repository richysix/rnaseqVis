#' Create UI components to select whether to cluster a matrix by rows and/or columns
#'
#' `clusterInput()` produces a checkbox group with check boxes for
#' clustering by genes (rows) and by samples (columns)
#' 
#' @param id namespace id for the UI components. Must match the id provided to the 
#' [clusterServer()] function.
#'
#' @returns a [shiny::checkboxGroupInput] object with two check boxes
#' 
#' @export
#'
#' @examples
#' 
#' clusterInput("rnaseqData")
#' 
clusterInput <- function(id) {
  tagList(
    checkboxGroupInput(
      NS(id, "clusterCheckGroup"),
      label = h4("Clustering"),
      choices = list("By Genes" = "genes", "By Samples" = "samples"),
      selected = c()
    )
  )
}

#' Server function to cluster a count matrix by rows and/or columns
#'
#' `clusterServer()` clusters the supplied counts wither by rows or columns
#' or both
#' 
#' @param id namespace id for the UI components. Must match the id provided to the 
#' [clusterInput()] function.
#' @param counts a reactive counts object. Should contain only numeric columns
#' @param gene_metadata a reactive object. Contains the metadata for the genes
#' present in the counts object.
#' @param debug Turn on debugging message statements
#'
#' @returns a [shiny::reactive()] object which is the clustered counts
#' 
#' @export
#'
#' @examples
#' 
#' clusterServer("rnaseqData", counts = reactive(rnaseqtools::counts[1:10,1:5]))
#' 
clusterServer <- function(id, counts = NULL, gene_metadata = NULL, debug = FALSE) {
  stopifnot(is.reactive(counts))
  stopifnot(is.reactive(gene_metadata))
  moduleServer(id, function(input, output, session) {
    counts_clustered <- reactive({
      req(counts())
      req(gene_metadata())
      
      if (debug) {
        print(counts())
        print(gene_metadata())
      }
      counts_m <- as.matrix(counts())
      if (debug) {
        print(counts_m[1:5,1:5])
        print(dim(counts_m))
        print(length(gene_metadata()$GeneID))
      }
      rownames(counts_m) <- gene_metadata()$GeneID
      
      possible_selections <- c("genes", "samples")
      index <- seq_along(possible_selections)[ possible_selections %in% input$clusterCheckGroup ] |> 
        sum()
      cluster_type <- c("none", "genes", "samples", "both")
      switch(
        cluster_type[index + 1],
        "none" = counts_m,
        "genes" = cluster_genes(counts_m),
        "samples" = cluster_samples(counts_m),
        "both" = cluster_both(counts_m)
      )
    })
    
    metadata_clustered <- reactive({
      req(counts_clustered())
      req(gene_metadata())
      # reorder metadata
      new_order <- sapply(rownames(counts_clustered()), 
                          function(x){ which(gene_metadata()$GeneID == x) }, 
                          USE.NAMES = FALSE)
      gene_metadata()[ new_order, ]
    })
      
    return(
      list(
        "counts" = counts_clustered,
        "gene_metadata" = metadata_clustered
      )
    )
  })
}

#' Cluster a gene count matrix
#'
#' `cluster_genes` clusters the matrix by rows (genes)
#' `cluster_samples` clusters the matrix by column (samples)
#' `cluster_both` clusters by both row and column
#' 
#' @param counts_m matrix of counts
#'
#' @return matrix 
#' @export
#'
#' @examples
#' counts_m <- matrix(rnorm(25), ncol = 5, 
#'   dimnames = list(paste("Gene", 1:5, sep = "-"), 
#'                   paste("sample", 1:5, sep = "-")))
#' 
#'  cluster_genes(counts_m)
#'  
#'  cluster_samples(counts_m)
#'  
#'  cluster_both(counts_m)
cluster_genes <- function(counts_m) {
  zeroVar <- RemoveZeroVariance(counts_m, rows = TRUE, cols = FALSE)
  miscr::cluster_matrix(zeroVar$matrix, by_row = TRUE, by_col = FALSE)
}

#' @rdname cluster_genes
cluster_samples <- function(counts_m) {
  zeroVar <- RemoveZeroVariance(counts_m, rows = FALSE, cols = TRUE)
  miscr::cluster_matrix(zeroVar$matrix, by_row = FALSE, by_col = TRUE)
}

#' @rdname cluster_genes
cluster_both <- function(counts_m) {
  zeroVar <- RemoveZeroVariance(counts_m, rows = TRUE, cols = TRUE)
  miscr::cluster_matrix(zeroVar$matrix, by_row = TRUE, by_col = TRUE)
}

#' RemoveZeroVariance
#'
#' \code{RemoveZeroVariance} removes rows/columns from a matrix that have zero variance
#'
#'    
#'    
#' @param x - matrix
#' @param rows - logical, default TRUE
#' @param cols - logical, default FALSE
#'
#' @return a list containing the following elements
#' 
#'    matrix - the supplied matrix with the rows/columns removed
#'    
#'    rowsKept - a character vector with the names of the rows that were not removed
#'      NULL if rows = FALSE
#'    
#'    rowsRemoved - a character vector with the names of the rows that were removed
#'      NULL if rows = FALSE
#'      
#'    colsKept - a character vector with the names of the columns that were not removed
#'      NULL if cols = FALSE
#'      
#'    colsRemoved - a character vector with the names of the columns that were removed
#'      NULL if cols = FALSE
#'
#' @examples
#' countMatrix <- matrix( sample(1:100, 100), ncol = 10)
#' countMatrix[ 3, ] <- rep(0,10)
#' countMatrix[ 6, ] <- rep(0,10)
#' RemoveZeroVariance( countMatrix, rows = TRUE, cols = TRUE )
#'
#' @export
#'
RemoveZeroVariance <- function( x, rows = TRUE, cols = FALSE ){
  # check matrix has row and col names
  if( rows & is.null(rownames(x)) ){
    stop("The supplied matrix does not have row names!")
  } else if( cols & is.null(colnames(x)) ){
    stop("The supplied matrix does not have column names!")
  }
  rowsKept <- NULL
  rowsRemoved <- NULL
  colsKept <- NULL
  colsRemoved <- NULL
  if( rows ){
    zeroVarRows <- genefilter::rowSds(x) == 0
    if (sum(zeroVarRows) > 0) {
      rowsKept <- rownames(x)[ !zeroVarRows ]
      rowsRemoved <- rownames(x)[ zeroVarRows ]
      x <- x[ !zeroVarRows, ]
    }
  } else{
    rowsKept <- rownames(x)
  }
  if( cols ){
    zeroVarCols <- genefilter::rowSds(t(x)) == 0
    if (sum(zeroVarCols) > 0) {
      colsKept <- colnames(x)[ !zeroVarCols ]
      colsRemoved <- colnames(x)[ zeroVarCols ]
      x <- x[ , !zeroVarCols ]
    }
  } else{
    colsKept <- colnames(x)
  }
  return(
    list(
      matrix = x,
      rowsKept = rowsKept,
      rowsRemoved = rowsRemoved,
      colsKept = colsKept,
      colsRemoved = colsRemoved
    )
  )
}

#' A test shiny app for the cluster module
#'
#' `clusterApp()` creates a small test app for testing the [clusterInput()] and
#' [clusterServer()] functions. It uses a subset of the package dataset `counts`,
#' clusters the counts based on the selected check boxes and displays them
#' as a table.
#' 
#' @return a [shiny::shinyApp()] object.
#'
#' @examples
#' 
#' clusterApp()
#' 
clusterApp <- function() {
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        clusterInput("cluster")
      ),
      mainPanel(
        verbatimTextOutput("rownames"),
        tableOutput("original_counts"),
        tableOutput("clustered_counts"),
        tableOutput("clustered_metadata")
      )
    )
  )
  
  server <- function(input, output, session) {
    test_counts <- rnaseqtools::counts[1:10, 1:5]
    counts_clustered <- clusterServer("cluster", counts = reactive(test_counts),
                                      gene_metadata = reactive(rnaseqtools::gene_metadata[1:10,]))
    output$rownames <- renderText(rownames(counts_clustered$counts()))
    output$original_counts <- renderTable(test_counts)
    output$clustered_counts <- renderTable(counts_clustered$counts())
    output$clustered_metadata <- renderTable(counts_clustered$metadata())
  }
  shinyApp(ui, server)
}
