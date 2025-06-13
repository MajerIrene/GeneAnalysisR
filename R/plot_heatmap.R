#' Plot a heatmap of gene expression grouped by clusters
#'
#' This function plots a heatmap of gene expression data, ordered by the clustering result.
#'
#' @param data A numeric matrix or data frame with genes as rows and samples as columns.
#' @param clustering_res The result from the `clustering()` function.
#' @param scale_data Logical, whether to scale rows before plotting (default TRUE).
#' @return NULL. The heatmap is plotted.
#'
#' @importFrom pheatmap pheatmap
#' @export
#' @examples
#' \dontrun{
#' hierarchical <- clustering(data, k = 10, method = "hierarchical", distance_method = 'euclidean')
#' plot_heatmap(normalized, hierarchical)
#' }
#'
plot_heatmap <- function(data, clustering_res, scale_data = TRUE) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Please install the 'pheatmap' package.")
  }

  if (!is.matrix(data)) data <- as.matrix(data)

  # Check if rownames exist and match cluster names
  if (is.null(rownames(data))) {
    warning("Data has no rownames; assigning generic rownames for plotting.")
    rownames(data) <- paste0("Gene", seq_len(nrow(data)))
  }

  if (is.null(names(clustering_res$clusters))) {
    warning("Cluster vector has no names; assuming order matches data rows.")
  } else if (!all(names(clustering_res$clusters) %in% rownames(data))) {
    warning("Not all cluster names found in data rownames; ordering may be incorrect.")
  }

  # Reorder genes by clusters
  ordered_idx <- order(clustering_res$clusters)
  ordered_data <- data[ordered_idx, , drop = FALSE]

  # Optional: scale rows (genes)
  if (scale_data) {
    ordered_data <- t(scale(t(ordered_data)))
  }

  annotation_row <- data.frame(Cluster = factor(clustering_res$clusters[ordered_idx]))
  rownames(annotation_row) <- rownames(ordered_data)

  pheatmap::pheatmap(
    ordered_data,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    annotation_row = annotation_row,
    main = "Heatmap of Clustered Gene Expression"
  )
}
