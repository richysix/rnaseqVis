#' Create UI components to upload sample and count files
#'
#' `uploadRNASeqInput()` produces the buttons for uploading sample and count data 
#' files. It also provides a checkbox to use the package test data. If sample and
#' count files are subsequently uploaded, this is used instead of the test data.
#' 
#' @param id namespace id for the UI components. Must match the id provided to the 
#' [uploadRNASeqServer()] function.
#'
#' @returns a [htmltools::tagList()] containing two [shiny::fileInput()] controls and 
#' a [shiny::checkboxInput()].
#' 
#' @export
#'
#' @examples
#' 
#' uploadRNASeqInput("rnaseqData")
#' 
uploadRNASeqInput <- function(id) {
  tagList(
    fileInput(NS(id, "sampleFile"), "Sample File"),
    fileInput(NS(id, "countFile"), "Count File"),
    checkboxInput(NS(id, "testdata"), 'Use test data', value = FALSE, width = NULL)
  )
}

#' Create Output to display uploaded sample and count data
#'
#' `uploadRNASeqOutput()` produces three table outputs for displaying 
#' the data returned by the `uploadRNASeqServer()` function. The three
#' outputs are labelled samples, counts and metadata. There are also two 
#' [shinyBS::bsAlert()] anchor points to alert the user to missing data in 
#' either the sample or count file
#' 
#' @param id namespace id for the UI components. Must match the id provided to the 
#' [uploadRNASeqServer()] function.
#'
#' @returns a [htmltools::tagList()] containing three [shiny::tableOutput()]s 
#' and two [shinyBS::bsAlert()] anchor points.
#' 
#' @export
#'
#' @examples
#' 
#' uploadRNASeqOutput("rnaseqData")
#' 
uploadRNASeqOutput <- function(id) {
  tagList(
    h3("Sample Data:"),
    shinyBS::bsAlert(NS(id, "countsInputAlert")),
    tableOutput("samples"),
    h3("Count Data:"),
    shinyBS::bsAlert(NS(id, "sampleInputAlert")),
    tableOutput("counts"),
    h3("Gene Metadata:"),
    tableOutput("metadata")
  )
}

#' Server function to upload sample and count data files
#'
#' `uploadRNASeqServer()` implements uploading a sample file and a count data
#' file. It also handles using the package test data.
#' 
#' @param id namespace id for the UI components. Must match the id provided to the 
#' [uploadRNASeqInput()] function.
#'
#' @returns a list containing two [shiny::reactive()] objects
#' * sampleInfo a data.frame of sample metadata
#' * counts a data.frame of RNAseq count data
#' 
#' @export
#'
#' @examples
#' 
#' uploadRNASeqServer("rnaseqData")
#' 
uploadRNASeqServer <- function(id, debug = FALSE) {
  moduleServer(id, function(input, output, session) {
    # set up observer to untick test data checkbox if a sample/count
    # file is uploaded
    observe({
      updateCheckboxInput(session, "testdata", value = FALSE)
    }) |>
      bindEvent(input$sampleFile, input$countFile)
    
    # return sample file path depending on whether the test data checkbox is checked
    sample_file <- reactive({
      if (input$testdata) {
        file_path <- system.file("extdata", "zfs-rnaseq-sampleInfo.tsv", package = "rnaseqVis")
        if(debug) message(paste0("Sample Info path = ", file_path))
        return(file_path)
      } else {
        if(debug) message(paste0("Sample Info path = ", input$sampleFile$name))
        return(input$sampleFile$datapath)
      }
    })
    
    # load sample data from file
    init_sample_info <- reactive({
      req(sample_file())
      rnaseqtools::load_rnaseq_samples(sample_file())
    })
    
    # return counts file path depending on whether the test data checkbox is checked
    counts_file <- reactive({
      if (input$testdata) {
        file_path <- system.file("extdata", "counts.shield-subset.tsv", package = "rnaseqVis")
        if(debug) message(paste0("Count Data path = ", file_path))
        return(file_path)
      } else {
        if(debug) message(paste0("Count Data path = ", input$countFile$name))
        return(input$countFile$datapath)
      }
    })
    # load rnaseq data from file
    init_rnaseq_data <- reactive({
      req(counts_file())
      counts <- rnaseqtools::load_rnaseq_data(counts_file())
      if (debug) {
        message("Loaded data file:")
        print(head(counts))
      }
      return(counts)
    })
    
    # This takes the samples and counts data frames and checks whether 
    # there are any samples in one that are not in the other
    # If there are, the extra samples are removed and a warning alert is created
    # The returned list contains sample and rnaseq_data data frames with the
    # same sample in each. This ensures that if we need to use DESeq2 to 
    # calculate normalised counts the sample data will match
    all_data <- reactive({
      req(init_sample_info())
      req(init_rnaseq_data())
      
      sample_info <- init_sample_info()
      rnaseq_data <- init_rnaseq_data()
      if (debug) {
        message("All data:")
        print(head(sample_info))
        print(head(rnaseq_data))
      }
      counts <- rnaseqtools::get_counts(rnaseq_data, normalised = TRUE)
      if (is.null(counts)) {
        counts <- rnaseqtools::get_counts(rnaseq_data)
      }
      
      # close any open alerts
      shinyBS::closeAlert(session, "sampleIdsAlert")
      shinyBS::closeAlert(session, "countSampleIdsAlert")
      sample_subset <- tryCatch(
        { rnaseqtools::check_samples_match_counts(counts, sample_info)
          return(sample_info) },
        warning = function(w) {
          if (any(grepl("missing_from.*_counts", class(w)))) {
            # subset sample info to samples in counts
            available_samples <- intersect(sample_info$sample, colnames(counts))
            if (length(available_samples) == 0){
              # create an alert and return NULL
              msg <- paste("<b>None</b> of the samples in the samples file match any of those in the counts file.",
                           sep = "<br>")
              # generate alert
              shinyBS::createAlert(session, anchorId = NS(id, "sampleInputAlert"),
                                   alertId = "sampleIdsAlert", title = "Sample IDs missing from counts",
                                   content = msg, append = FALSE, style = "danger"
              )
              return(NULL)
            }
            samples <- sample_info[ sample_info$sample %in% available_samples, ]
            # parse message
            if ("missing_from_both_samples_and_counts" %in% class(w)) {
              # remove second half of message
              msg <- sub("\n.*$", "", w$message)
            } else {
              msg <- w$message
            }
            msg <- paste(msg, "These samples have been removed from the sample information",
                         "If you want these samples included, they must also be present in the counts file",
                         sep = "<br>")
            # generate alert
            shinyBS::createAlert(session, anchorId = NS(id, "sampleInputAlert"),
                                 alertId = "sampleIdsAlert", title = "Sample IDs missing from counts",
                                 content = msg, append = FALSE, style = "warning"
            )
          } else {
            samples <- sample_info
          }
          if (debug) {
            message("All data: warning from check_samples_match_counts")
            message(w$message)
          }
          return(samples)
        }
      )
      if (debug) {
        message("All data: subset samples")
        print(sample_subset)
      }
      rnaseq_data_subset <- tryCatch(
        { rnaseqtools::check_samples_match_counts(counts, sample_info)
          return(rnaseq_data) },
        warning = function(w) {
          if (any(grepl("missing_from.*_samples", class(w)))) {
            # subset rnaseq_data to samples
            rnaseq_subset <- rnaseqtools::subset_to_samples(rnaseq_data, sample_info)
            count_samples <- colnames(rnaseq_data) |> 
              (function(x){ sub(" count", "", x) })() |> 
              (function(x){ sub(" normalised", "", x) })()
            matching_samples <- intersect(sample_info$sample, count_samples)
            if (length(matching_samples) == 0){
              # create an alert and return NULL
              msg <- paste("<b>None</b> of the samples in the counts file match any of those in the samples file.",
                           sep = "<br>")
              # generate alert
              shinyBS::createAlert(session, anchorId = NS(id, "countsInputAlert"),
                                   alertId = "countSampleIdsAlert", title = "Sample IDs missing from samples file",
                                   content = msg, append = FALSE, style = "danger"
              )
              return(NULL)
            }
            # parse message
            if ("missing_from_both_samples_and_counts" %in% class(w)) {
              # remove first half of message
              msg <- sub("^.*\n", "", w$message) |> 
                (function(x){ sub(" Only samples in both were returned", "", x) })()
            } else {
              msg <- w$message
            }
            msg <- paste(msg, "These samples have been removed from the count data",
                         "If you want these samples included, they must be present in the samples file",
                         sep = "<br>")
            # generate alert
            shinyBS::createAlert(session, anchorId = NS(id, "countsInputAlert"),
                                 alertId = "countSampleIdsAlert", title = "Sample IDs missing from samples file",
                                 content = msg, append = FALSE, style = "warning"
            )
          } else {
            rnaseq_subset <- rnaseq_data
          }
          if (debug) message(w$message)
          return(rnaseq_subset)
        }
      )
      if (debug) {
        message("All data: subset RNAseq data")
        print(rnaseq_data_subset)
      }
      list(
        "sample_info" = sample_subset,
        "rnaseq_data" = rnaseq_data_subset
      )
    })
    
    # unpack the all_data reactive list
    rnaseq_data <- reactive({
      req(all_data())
      all_data()$rnaseq_data
    })
    sample_info <- reactive({
      req(all_data())
      all_data()$sample_info
    })
    
    # get normalised counts, either from the rnaseq data object or
    # by using DESeq2 to calculate them
    # This may need a progress bar at some point
    count_data <- reactive({
      req(rnaseq_data())
      req(sample_info())
      norm_counts <- rnaseqtools::get_counts(rnaseq_data(), normalised = TRUE)
      # if file does not have normalised counts, get raw counts and normalise
      if (is.null(norm_counts)) {
        norm_data <- rnaseqtools::normalise_counts(rnaseq_data(), sample_info())
        if (debug) {
          message("count_data reactive: normalise_counts")
          print(head(norm_data))
        }
        norm_counts <- rnaseqtools::get_counts(norm_data, normalised = TRUE)
      }
      return(norm_counts)
    })
    
    # extract gene metadata
    gene_metadata <- reactive({
      req(rnaseq_data())
      rnaseqtools::get_gene_metadata(rnaseq_data())
    })
    # return the sample info and count data
    list(
      sample_info = sample_info,
      gene_metadata = gene_metadata,
      counts = count_data
    )
  })

}

#' A test shiny app for the uploadRNASeq module
#'
#' `uploadRNASeqApp()` creates a small test app for testing the [uploadRNASeqInput()] and
#' [uploadRNASeqServer()] functions. A subset of the returned sample and count data.frames 
#' are displayed in two [shiny::tableOutput()]s
#' 
#' @return a [shiny::shinyApp()] object
#'
#' @examples
#' uploadRNASeqApp()
uploadRNASeqApp <- function(debug = FALSE) {
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
        uploadRNASeqInput("rnaseqData")
      ),
      mainPanel(
        uploadRNASeqOutput("rnaseqData"),
      )
    )
  )
  server <- function(input, output, session) {
    data_list <- uploadRNASeqServer("rnaseqData", debug)
    output$samples <- renderTable(data_list$sample_info()[1:5,1:5])
    output$counts <- renderTable(data_list$counts()[1:5,1:10])
    output$metadata <- renderTable(data_list$gene_metadata()[1:5,])
  }
  shinyApp(ui, server)
}

