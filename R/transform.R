#' Create UI components to select a transformation function for count data
#'
#' `transformInput()` produces a set of radio buttons for selecting a function
#' to use to transform the count data
#' 
#' @param id namespace id for the UI components. Must match the id provided to the 
#' [transformServer()] function.
#'
#' @returns a [shiny::radioButtons()] object with four choices
#' 
#' @export
#'
#' @examples
#' 
#' transformInput("rnaseqData")
#' 
transformInput <- function(id) {
  radioButtons(
    NS(id, "transform_func"),
    label = h4("Transform Counts"),
    choices = list(
      "Raw" = "raw",
      "Max Scaled" = "max",
      "log 10" = "log",
      "Mean Centred and Scaled" = "zscore"
    ),
    selected = "raw"
  )
}

#' Server function to transform a count matrix with a selected function
#'
#' `transformServer()` transforms the supplied counts with the appropriate function
#' based on the radio buttons selection
#' 
#' @param id namespace id for the UI components. Must match the id provided to the 
#' [transformInput()] function.
#' @param counts a reactive counts object. Should contain only numeric columns
#'
#' @returns a list of two [shiny::reactive()] objects
#' * counts - the transformed counts
#' * transform - the name of the selected transform
#' 
#' @export
#'
#' @examples
#' 
#' transformServer("rnaseqData", counts = reactive(rnaseqtools::counts[1:10,1:5]))
#' 
transformServer <- function(id, counts = NULL) {
  stopifnot(is.reactive(counts))
  moduleServer(id, function(input, output, session) {
    counts_transformed <- reactive({
      req(input$transform_func)
      if (is.null(counts())) {
        return(NULL)
      } else {
        return(transform_counts(counts(), input$transform_func))
      }
    })

    return(
      list(
        "counts" = counts_transformed,
        "transform" = reactive(input$transform_func)
      )
    )
  })
}

#' Transform counts with the supplied function name
#'
#' @param counts tibble with only numeric columns
#' @param transform_func character string to select the transform function
#'
#' @return a [tibble::tibble] conatining the transformed counts
#' @export
#'
#' @examples
#' 
#' transform_counts(rnaseqtools::counts, "raw")
#' 
#' transform_counts(rnaseqtools::counts, "zscore")
#' 
transform_counts <- function(counts, transform_func) {
  transformed_counts <- switch(
    transform_func,
    "raw" = counts,
    "max" = max_scale(counts),
    "log" = log10(counts + 1) |> tibble::as_tibble(),
    "zscore" = t(scale(t(counts))) |> tibble::as_tibble()
  )
  return(transformed_counts)
}

#' Scale count according to the maximum value
#'
#' @param counts tibble of counts
#'
#' @return a tibble of max sclaed counts
#' @export
#'
#' @examples
#' 
#' max_scale(rnaseqtools::counts)
#' 
max_scale <- function(counts) {
  geneMaxCounts <- apply(counts, 1, max)
  # scale operates on the column so need to transpose, scale and then transpose back
  counts_max_scaled <-
    t(scale(t(counts), scale = geneMaxCounts, center = FALSE))
  # genes with a max of zero get converted to NAs
  # reset to zeros
  counts_max_scaled[ geneMaxCounts == 0, ] <- 
    matrix( rep(0, sum(geneMaxCounts == 0)*ncol(counts_max_scaled) ), ncol = ncol(counts_max_scaled) )
  # remove attributes introduced by scale
  return(tibble::as_tibble(counts_max_scaled[,]))
}

#' A test shiny app for the transform module
#'
#' `transformApp()` creates a small test app for testing the [transformInput()] and
#' [transformServer()] functions. It uses a subset of the package dataset `counts`,
#' transforms the counts based on the selected radio buttons and displays them
#' as a table.
#' 
#' @return a [shiny::shinyApp()] object.
#'
#' @examples
#' 
#' transformApp()
#' 
transformApp <- function() {
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        transformInput("transform")
      ),
      mainPanel(
        verbatimTextOutput('func'),
        tableOutput("transformed_counts")
      )
    )
  )
  
  server <- function(input, output, session) {
    test_counts <- rnaseqtools::counts[1:10, 1:5]
    counts_transformed <- transformServer("transform", counts = reactive(test_counts))
    output$func <- renderText(counts_transformed$transform())
    output$transformed_counts <- renderTable(counts_transformed$counts())
  }
  shinyApp(ui, server)
} 
