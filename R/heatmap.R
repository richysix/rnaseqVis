#' Create UI components to display heatmap of a matrix of values
#'
#' `heatmapOutput()` produces a plotOutput space for a plot
#' 
#' @param id namespace id for the UI components. Must match the id provided to the 
#' [heatmapServer()] function.
#'
#' @returns a [htmltools::tagList()] containing a [shiny::plotOutput] object
#' 
#' @export
#'
#' @examples
#' 
#' heatmapOutput("rnaseqData")
#' 
heatmapOutput <- function(id) {
  tagList(
    plotOutput(NS(id, "heatmap_plot"))
  )
}

#' Server function to create a heatmap of a count matrix
#'
#' `heatmapServer()` creates a [ggplot2] plot object using [biovisr::matrix_heatmap]
#' from the supplied count matrix, sample info and gene metadata
#' 
#' @param id namespace id for the UI components. Must match the id provided to the 
#' [heatmapInput()] function.
#' @param counts a reactive counts object. Should contain only numeric columns
#' @param sample_info a reactive object. Represents the samples and associated
#' metadata
#' @param gene_metadata a reactive object. Contains the metadata for the genes
#' present in the counts object.
#' @param transform reactive that contains the contents of input$transform from
#' the transform module
#' @param debug Turn on debugging message statements
#'
#' @returns a [shiny::reactive()] object which is the heatmaped counts
#' 
#' @export
#'
#' @examples
#' 
#' heatmapServer("rnaseqData", counts = reactive(rnaseqtools::counts[1:10,1:5]))
#' 
heatmapServer <- function(id, counts = NULL, sample_info = NULL,
                          gene_metadata = NULL, transform = NULL,
                          debug = NULL) {
  stopifnot(is.reactive(counts))
  stopifnot(is.reactive(sample_info))
  stopifnot(is.reactive(gene_metadata))
  stopifnot(is.reactive(transform))

  moduleServer(id, function(input, output, session) {
    plot <- reactive({
      req(counts())
      req(sample_info())
      req(gene_metadata())
      req(transform())

      # show gene names if matrix is small enough
      counts <- counts()
      if (nrow(counts) <= 100) {
        gene_names <- get_gene_labels(gene_metadata())
      } else {
        gene_names <- FALSE
      }
      # show sample names if matrix is small enough
      sample_info <- sample_info()
      if (ncol(counts) <= 48) {
        sample_names <- get_sample_labels(sample_info(), colnames(counts))
      } else {
        sample_names <- FALSE
      }

      plot <- biovisr::matrix_heatmap(counts, xaxis_labels = sample_names, yaxis_labels = gene_names)
      if (transform() == 'zscore') {
        plot <- plot +
          ggplot2::scale_fill_distiller(type = 'div', palette = "RdBu")
      }

      return(plot)
    })

    output$heatmap_plot <- renderPlot(plot())
  })
}

#' Get labels for heatmap axis labels
#'
#' @param gene_metadata data.frame - Gene metadata
#' @param sample_info data.frame - Sample info
#'
#' @return vector of labels
#' @export
#'
#' @examples
#' 
#' get_gene_labels(gene_metadata)
#' 
#' get_sample_labels(sample_info)
#' 
get_gene_labels <- function(gene_metadata) {
  if ("Name" %in% names(gene_metadata)) {
    return(gene_metadata$Name)
  } else {
    return(gene_metadata$GeneID)
  }
}

#' @rdname get_gene_labels
get_sample_labels <- function(sample_info, sample_ids) {
  if ("sampleName" %in% names(sample_info)) {
    labels <- sample_info$sampleName
    names(labels) <- sample_info$sample
    return(labels[ sample_ids ])
  } else {
    return(sample_info$sample)
  }
}

#' A test shiny app for the heatmap module
#'
#' `heatmapApp()` creates a small test app for testing the [heatmapOutput()] and
#' [heatmapServer()] functions. It uses a subset of the package datasets `counts`,
#' `sampleInfo` and `gene_metadata` and create a [ggplot2::ggplot()] heatmap object.
#' It also has transform radio buttons to test using that information to change
#' the colour palette.
#' 
#' @return a [shiny::shinyApp()] object.
#'
#' @examples
#' 
#' heatmapApp()
#' 
heatmapApp <- function() {
  ui <- fluidPage(
    heatmapOutput('heatmap'),
    radioButtons(
      "transform_func",
      label = h4("Transform Counts"),
      choices = list(
        "Raw" = "raw",
        "Max Scaled" = "max",
        "log 10" = "log",
        "Mean Centred and Scaled" = "zscore"
      ),
      selected = "raw"
    )
  )

  server <- function(input, output, session) {
    heatmapServer("heatmap", counts = reactive(rnaseqtools::counts[1:10, 1:5]),
                  sample_info = reactive(rnaseqtools::sampleInfo[1:5,]),
                  gene_metadata = reactive(rnaseqtools::gene_metadata[1:10,]),
                  transform = reactive(input$transform_func),
                  debug = TRUE)
  }
  shinyApp(ui, server)
}
