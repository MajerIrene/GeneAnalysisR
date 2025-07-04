% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_dendrogram.R
\name{plot_dendrogram}
\alias{plot_dendrogram}
\title{Plot improved dendrogram from clustering result}
\usage{
plot_dendrogram(clustering_res)
}
\arguments{
\item{clustering_res}{A list returned by the `clustering` function,
with an element `model` of class `hclust` and `clusters` vector.}
}
\value{
Plots the dendrogram.
}
\description{
This function plots a hierarchical clustering dendrogram with improved
visualization: smaller, colored labels and horizontal orientation.
}
\examples{
\dontrun{
hierarchical <- clustering(data, k = 4, method = "hierarchical")
plot_dendrogram(hierarchical)
}
}
