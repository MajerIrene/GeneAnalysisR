#' Plot improved dendrogram from clustering result
#'
#' This function plots a hierarchical clustering dendrogram with improved
#' visualization: smaller, colored labels and horizontal orientation.
#'
#' @param clustering_res A list returned by the `clustering` function,
#'        with an element `model` of class `hclust` and `clusters` vector.
#' @return Plots the dendrogram.
#'
#' @importFrom stats as.dendrogram
#' @importFrom dendextend color_branches color_labels set rect.dendrogram
#' @export
#' @examples
#' \dontrun{
#' hierarchical <- clustering(data, k = 4, method = "hierarchical")
#' plot_dendrogram(hierarchical)
#' }

plot_dendrogram <- function(clustering_res) {
  if (!is.list(clustering_res) || is.null(clustering_res$model) || is.null(clustering_res$clusters)) {
    stop("Input must be a list with 'model' and 'clusters'.")
  }
  if (!inherits(clustering_res$model, "hclust")) {
    stop("Model is not hierarchical.")
  }

  dend <- as.dendrogram(clustering_res$model)
  k <- length(unique(clustering_res$clusters))

  dend <- color_branches(dend, k = k)
  dend <- set(dend, "labels", NA)  # Labels are not readable so I won't show them

  old_mar <- par("mar")
  on.exit(par(mar = old_mar))
  par(mar = c(5, 4, 4, 2))

  plot(dend, main = "Hierarchical Clustering Dendrogram", horiz = TRUE)
  rect.dendrogram(dend, k = k, border = 2:(k + 1), horiz = TRUE)
}
