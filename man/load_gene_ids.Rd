% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SubsetByGeneIDs.R
\name{load_gene_ids}
\alias{load_gene_ids}
\title{Load a file of gene ids}
\usage{
load_gene_ids(gene_ids_file)
}
\arguments{
\item{gene_ids_file}{path to file of gene ids.}
}
\value{
a vector of the first column of the file.
If \code{gene_ids_file} is NULL, then it returns NULL
}
\description{
\code{load_gene_ids()} takes a file name and returns the first column of the
file.
}
\examples{
ids_file <- system.file("extdata", "gene-ids.txt", package = "rnaseqVis")
load_gene_ids(ids_file)

}
