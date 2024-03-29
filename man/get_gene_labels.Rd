% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap.R
\name{get_gene_labels}
\alias{get_gene_labels}
\alias{get_sample_labels}
\title{Get labels for heatmap axis labels}
\usage{
get_gene_labels(gene_metadata)

get_sample_labels(sample_info, sample_ids)
}
\arguments{
\item{gene_metadata}{data.frame - Gene metadata}

\item{sample_info}{data.frame - Sample info}
}
\value{
vector of labels
}
\description{
Get labels for heatmap axis labels
}
\examples{

get_gene_labels(gene_metadata)

get_sample_labels(sample_info)

}
