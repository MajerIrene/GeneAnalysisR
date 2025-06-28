#' Plot improved dendrogram from clustering result
#'
#' This function plots a hierarchical clustering dendrogram with improved
#' visualization: smaller, colored labels and horizontal orientation to improve readability.
#'
#' @param clustering_res A list returned by the `clustering` function,
#'        with an element `model` of class `hclust` and `clusters` vector.
#' @return NULL (plots the dendrogram).
#'
#' @importFrom stats as.dendrogram
#' @importFrom dendextend color_branches color_labels set rect.dendrogram
#' @export
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(100), nrow = 50)
#' rownames(data) <- paste0("Gene", 1:50)
#' res <- clustering(data, k = 4, method = "hierarchical")
#' plot_dendrogram_improved(res)
#' }

plot_dendrogram <- function(clustering_res) {
  if (!requireNamespace("dendextend", quietly = TRUE)) {
    stop("Please install the 'dendextend' package.")
  }

  if (!is.list(clustering_res) || is.null(clustering_res$model) || is.null(clustering_res$clusters)) {
    stop("Input must be a list with 'model' and 'clusters'.")
  }
  if (!inherits(clustering_res$model, "hclust")) {
    stop("Model is not hierarchical.")
  }

  dend <- as.dendrogram(clustering_res$model)
  k <- length(unique(clustering_res$clusters))

  dend <- color_branches(dend, k = k)
  dend <- set(dend, "labels", NA)  # rimuove etichette

  old_mar <- par("mar")
  on.exit(par(mar = old_mar))
  par(mar = c(5, 4, 4, 2))  # meno spazio per etichette

  plot(dend, main = "Hierarchical Clustering Dendrogram", horiz = TRUE)
  rect.dendrogram(dend, k = k, border = 2:(k + 1), horiz = TRUE)
}
