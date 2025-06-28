#' Plot a heatmap of gene expression grouped by clusters
#'
#' This function plots a heatmap of gene expression data, ordered by the clustering result.
#'
#' @param expr_data A numeric matrix of raw counts (genes x samples).
#' @param clustering_res The result from the `clustering()` function.
#' @param scale_data Logical, whether to scale rows before plotting. Default is TRUE.
#'
#' @return The heatmap is plotted.
#'
#' @importFrom pheatmap pheatmap
#' @export
#'
#' @examples
#' \dontrun{
#' hierarchical <- clustering(normalized, k = 10, method = "hierarchical", distance_method = 'euclidean')
#' plot_heatmap(normalized, hierarchical)
#' }
plot_heatmap <- function(expr_data, clustering_res, scale_data = TRUE) {

  if (!is.matrix(expr_data)) expr_data <- as.matrix(expr_data)
  # Adding check on input after test
  if (is.null(clustering_res) || !is.list(clustering_res)) {
    stop("clustering_res must be a list.")
  }

  if (is.null(clustering_res$model)) {
    stop("clustering_res$model is missing.")
  }

  if (is.null(clustering_res$clusters)) {
    stop("clustering_res$clusters is missing or NULL.")
  }

  if (!is.vector(clustering_res$clusters)) {
    stop("clustering_res$clusters must be a vector.")
  }

  if (length(clustering_res$clusters) != nrow(expr_data)) {
    stop("Length of clustering_res$clusters does not match number of rows in expr_data.")
  }

  # Check if rownames exist and match cluster names
  if (is.null(rownames(expr_data))) {
    warning("expr_data has no rownames; assigning generic rownames for plotting.")
    rownames(expr_data) <- paste0("Gene", seq_len(nrow(expr_data)))
  }

  if (is.null(names(clustering_res$clusters))) {
    warning("Cluster vector has no names; assuming order matches expr_data rows.")
  } else if (!all(names(clustering_res$clusters) %in% rownames(expr_data))) {
    warning("Not all cluster names found in expr_data rownames; ordering may be incorrect.")
  }

  # Reorder genes by clusters
  ordered_idx <- order(clustering_res$clusters)
  ordered_data <- expr_data[ordered_idx, , drop = FALSE]

  # Optionally scale rows
  if (scale_data) {
    ordered_data <- t(scale(t(ordered_data)))
  }

  annotation_row <- data.frame(Cluster = factor(clustering_res$clusters[ordered_idx]))
  rownames(annotation_row) <- rownames(ordered_data)

  pheatmap(
    ordered_data,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    annotation_row = annotation_row,
    main = "Heatmap of Clustered Gene Expression"
  )
}
