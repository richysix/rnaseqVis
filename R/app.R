library(shiny)
library(shinyjs)

rnaseqVisApp <- function(debug = FALSE, ...) {
  ui <- bslib::page_navbar(
    title = "rnaseqVis",
    shinyjs::useShinyjs(),
    theme = bslib::bs_theme(bootswatch = "cosmo"),
    tabPanel(
      "Input",
      sidebarLayout(
        sidebarPanel(
          shinyModules::uploadRNASeqInput("upload"),
          SubsetByGeneIDsInput("subset"),
          width = 3
        ),
        mainPanel(
          SubsetByGeneIDsOutput("subset"),
          shinyModules::uploadRNASeqOutput("upload"),
        )
      )
    ),
    tabPanel(
      "Heatmap",
      fluidRow(
        transformInput("transform"),
        clusterInput("cluster"),
      ),
      heatmapOutput("hmap")
    ),
    tabPanel(
      "Count Plot",
      sidebarLayout(
        sidebarPanel(
          shinyModules::countPlotInput('countPlot'),
          width = 3
        ),
        mainPanel(
          fluidRow(
            shinyModules::countPlotOutput('countPlot'),
          ),
          width = 9
        )
      )
    ),
    tabPanel("Help", includeMarkdown("README.md"))
  )
  server <- function(input, output, session) {
    # increase max upload size to 200 MB
    options(shiny.maxRequestSize = 200 * 1024 ^ 2)
    
    data_list <- shinyModules::uploadRNASeqServer("upload", debug = debug)

    counts_subset <- SubsetByGeneIDsServer(
      "subset", 
      "counts" = reactive(data_list$counts()),
      "gene_metadata" = reactive(data_list$gene_metadata()),
      debug
    )

    # Transform counts
    transformed_counts <- transformServer("transform", counts = reactive(counts_subset$counts_subset()))
    
    # Cluster counts
    clustered_counts <- clusterServer(
      "cluster", counts = reactive(transformed_counts$counts()), 
      gene_metadata = reactive(counts_subset$gene_metadata_subset()),
      debug = debug
    )
    
    heatmapServer(
      "hmap", counts = reactive(clustered_counts$counts()),
      sample_info = reactive(data_list$sample_info()),
      gene_metadata = reactive(clustered_counts$gene_metadata()),
      transform = reactive(transformed_counts$transform()),
      debug = debug
    )
    
    shinyModules::countPlotServer(
      "countPlot", counts = reactive(counts_subset$counts()), 
      sample_info = reactive(data_list$sample_info()),
      gene_metadata = reactive(counts_subset$gene_metadata_subset()),
      debug = debug
    )
    
    # Outputs
    output$samples <- renderTable(head(data_list$sample_info()))
    output$counts <- renderTable(clustered_counts$counts()[1:5,1:5])
    output$metadata <- renderTable(clustered_counts$gene_metadata()[1:5,])
  }
  shinyApp(ui, server, ...)
}
