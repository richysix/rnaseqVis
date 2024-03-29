% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.R
\name{RemoveZeroVariance}
\alias{RemoveZeroVariance}
\title{RemoveZeroVariance}
\usage{
RemoveZeroVariance(x, rows = TRUE, cols = FALSE)
}
\arguments{
\item{x}{\itemize{
\item matrix
}}

\item{rows}{\itemize{
\item logical, default TRUE
}}

\item{cols}{\itemize{
\item logical, default FALSE
}}
}
\value{
a list containing the following elements

matrix - the supplied matrix with the rows/columns removed

rowsKept - a character vector with the names of the rows that were not removed
NULL if rows = FALSE

rowsRemoved - a character vector with the names of the rows that were removed
NULL if rows = FALSE

colsKept - a character vector with the names of the columns that were not removed
NULL if cols = FALSE

colsRemoved - a character vector with the names of the columns that were removed
NULL if cols = FALSE
}
\description{
\code{RemoveZeroVariance} removes rows/columns from a matrix that have zero variance
}
\examples{
countMatrix <- matrix( sample(1:100, 100), ncol = 10)
countMatrix[ 3, ] <- rep(0,10)
countMatrix[ 6, ] <- rep(0,10)
RemoveZeroVariance( countMatrix, rows = TRUE, cols = TRUE )

}
