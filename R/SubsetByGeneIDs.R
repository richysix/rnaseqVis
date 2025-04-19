#' Create UI components to upload a file to use to subset the data
#' 
#' `SubsetByGeneIDsInput()` produces upload and reset buttons to provide a file of 
#' gene IDs to use to subset the data with. 
#' 
#' @param id namespace id for the UI components. Must match the id provided to the 
#' [SubsetByGeneIDsServer()] function.
#' 
#' @returns a [htmltools::tagList()] containing a [shiny::fileInput()] control and 
#' a [shiny::actionButton()].
#' 
#' @examples 
#' SubsetByGeneIDsInput("geneIds")
#' 
#' @export
SubsetByGeneIDsInput <- function(id) {
  tagList(
    fileInput(NS(id, "geneIdsFile"), "Upload file to subset genes with"),
    actionButton(NS(id, "subsetReset"), "Reset Subset Gene IDs"),
  )
}

#' Create UI components to upload a file to use to subset the data
#' 
#' `SubsetByGeneIDsOutput()` creates a [shiny::uiOutput()] to
#' allow alerting the user if any of the supplied gene IDs are not found in
#' the data
#' 
#' @param id namespace id for the UI components. Must match the id provided to the 
#' [SubsetByGeneIDsServer()] function.
#' 
#' @returns a [htmltools::tagList()] containing a [shiny::uiOutput()] anchor point
#' 
#' @examples 
#' SubsetByGeneIDsOutput("geneIds")
#' 
#' @export
SubsetByGeneIDsOutput <- function(id) {
  tagList(
    uiOutput(NS(id, "subsetAlert")),
  )
}

#' Server function to upload a file to use to subset the data
#' 
#' `SubsetByGeneIDsServer()` implements uploading a file of gene IDs to use to 
#' subset the data with. It also deals with resetting the gene IDs with the 
#' reset button.
#' 
#' @param id namespace id for the UI components. Must match the id provided to the 
#' [SubsetByGeneIDsInput()] function.
#' @param counts a reactive counts object. Should contain only numeric columns
#' @param gene_metadata a reactive object. Contains the metadata for the genes
#' present in the counts object.
#' @param debug Turn on debugging message statements
#' 
#' @returns a [shiny::reactive()] object which is a vector of gene IDs if a file has been uploaded.
#' If a file hasn't been uploaded or the reset button has been clicked it returns NULL
#' 
#' @examples 
#' SubsetByGeneIDsServer("geneIds")
#' 
#' # turn on debugging
#' SubsetByGeneIDsServer("geneIds", debug = TRUE)
#' 
#' @export
#' 
SubsetByGeneIDsServer <- function(id, counts = NULL, gene_metadata = NULL, debug = FALSE) {
  stopifnot(is.reactive(counts))
  stopifnot(is.reactive(gene_metadata))
  stopifnot(is.logical(debug))
  moduleServer(id, function(input, output, session) {
    upload_state <- reactiveValues(
      state = NULL
    )
    observe({
      if (debug) message("Gene IDs file uploaded")
      upload_state$state <- 'uploaded'
    }, label = "Upload state") |> bindEvent(input$geneIdsFile)
    observe({
      if (debug) message("Gene IDs file reset")
      upload_state$state <- 'reset'
    }, label = "Reset state") |> bindEvent(input$subsetReset)
    
    geneIdsFile <- reactive({
      if (is.null(input$geneIdsFile) | is.null(upload_state$state)) {
        return(NULL)
      }
      if (upload_state$state == "uploaded") {
        return(input$geneIdsFile$datapath)
      } else {
        return(NULL)
      }
    })
    gene_ids <- reactive({
      load_gene_ids(geneIdsFile())
    })
    
    subset <- reactive({
      req(counts())
      req(gene_metadata())
      
      if (debug) {
        print(head(counts()))
        print(head(gene_metadata()))
        print(gene_ids())
      }
      if (is.null(gene_ids())) {
        counts_list <- list(
            "counts_subset" = counts(),
            "gene_metadata_subset" = gene_metadata()
        )
      } else {
        # close any open alerts
        shinyjs::runjs(glue::glue('$("#{id}-all_genes_missing button").click()'))
        shinyjs::runjs(glue::glue('$("#{id}-genes_missing button").click()'))
        # find any genes in gene_ids that don't exist in counts
        # first check whether none of them exist
        if (all(!(gene_ids() %in% gene_metadata()$GeneID))) {
          missing_genes <- setdiff(gene_ids(), gene_metadata()$GeneID)
          output$subsetAlert <- renderUI(
            shinyWidgets::alert(
              tags$h4("Gene IDs missing from counts"),
              tags$b("None"), "of the supplied gene IDs were found in the data:",
              tags$br(), paste0(missing_genes, collapse = ", "),
              tags$br(), "The original data has been returned",
              status = "danger",
              dismissible = TRUE
            )
          )
          selected_rows <- rep(TRUE, length(gene_metadata()$GeneID))
        } else if (any(!(gene_ids() %in% gene_metadata()$GeneID))) {
          # create alert
          missing_genes <- setdiff(gene_ids(), gene_metadata()$GeneID)
          output$subsetAlert <- renderUI(
            shinyWidgets::alert(
              tags$h4("Gene IDs missing from counts"),
              "The following gene IDs where not found in the data to subset:",
              tags$br(), paste0(missing_genes, collapse = ", "),
              status = "danger",
              dismissible = TRUE
            )
          )
          selected_rows <- gene_metadata()$GeneID %in% gene_ids()
        } else {
          selected_rows <- gene_metadata()$GeneID %in% gene_ids()
        }
        
        counts_list <- list(
          "counts_subset" = counts()[ selected_rows, ],
          "gene_metadata_subset" = gene_metadata()[ selected_rows, ]
        )
      }
      return(counts_list)
    })
    
    list(
      "counts_subset" = reactive(subset()$counts_subset),
      "gene_metadata_subset" = reactive(subset()$gene_metadata_subset)
    )
  })
}

#' Load a file of gene ids
#'
#' `load_gene_ids()` takes a file name and returns the first column of the
#' file.
#' 
#' @param gene_ids_file path to file of gene ids.
#'
#' @returns a vector of the first column of the file.
#' If `gene_ids_file` is NULL, then it returns NULL
#' 
#' @examples
#' ids_file <- system.file("extdata", "gene-ids.txt", package = "rnaseqVis")
#' load_gene_ids(ids_file)
#' 
load_gene_ids <- function(gene_ids_file) {
  if(is.null(gene_ids_file)){
    return(NULL)
  } else {
    data <- readr::read_tsv(gene_ids_file, show_col_types = FALSE,
                            col_names = FALSE)
    return(data[[1]])
  }
}

#' A test shiny app for the SubsetByGeneIDs module
#'
#' `SubsetByGeneIDsApp()` creates a small test app for testing the [SubsetByGeneIDsInput()] and
#' [SubsetByGeneIDsServer()] functions. A subset of the returned gene IDs and count data.frames 
#' are displayed in two [shiny::tableOutput()]s
#' 
#' @param debug logical for turning on debugging messages. default: FALSE
#' 
#' @return a [shiny::shinyApp()] object.
#'
#' @examples
#' SubsetByGeneIDsApp()
#' 
#' SubsetByGeneIDsApp(debug = TRUE)
#' 
SubsetByGeneIDsApp <- function(debug = FALSE) {
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        SubsetByGeneIDsInput("geneIds")
      ),
      mainPanel(
        SubsetByGeneIDsOutput("geneIds"),
        h3("Count Data:"),
        tableOutput('counts'),
        h3("Gene Metadata:"),
        tableOutput('metadata')
      )
    )
  )
  
  server <- function(input, output, session) {
    data_list <- SubsetByGeneIDsServer("geneIds", 
                                       "counts" = reactive(rnaseqtools::counts),
                                       "gene_metadata" = reactive(rnaseqtools::gene_metadata),
                                       "debug" = debug)
    output$counts <- renderTable(data_list$counts_subset()[1:5,1:6])
    output$metadata <- renderTable(data_list$gene_metadata_subset()[1:5,])
  }
  shinyApp(ui, server)
} 
