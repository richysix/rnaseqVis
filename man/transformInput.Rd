% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform.R
\name{transformInput}
\alias{transformInput}
\title{Create UI components to select a transformation function for count data}
\usage{
transformInput(id)
}
\arguments{
\item{id}{namespace id for the UI components. Must match the id provided to the
\code{\link[=transformServer]{transformServer()}} function.}
}
\value{
a \code{\link[shiny:radioButtons]{shiny::radioButtons()}} object with four choices
}
\description{
\code{transformInput()} produces a set of radio buttons for selecting a function
to use to transform the count data
}
\examples{

transformInput("rnaseqData")

}
