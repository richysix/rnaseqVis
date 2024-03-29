% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap.R
\name{heatmapOutput}
\alias{heatmapOutput}
\title{Create UI components to display heatmap of a matrix of values}
\usage{
heatmapOutput(id)
}
\arguments{
\item{id}{namespace id for the UI components. Must match the id provided to the
\code{\link[=heatmapServer]{heatmapServer()}} function.}
}
\value{
a \code{\link[htmltools:tagList]{htmltools::tagList()}} containing a \link[shiny:plotOutput]{shiny::plotOutput} object
}
\description{
\code{heatmapOutput()} produces a plotOutput space for a plot
}
\examples{

heatmapOutput("rnaseqData")

}
