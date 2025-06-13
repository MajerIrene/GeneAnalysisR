#' Compute log2 Fold Change for Genes in Each Cluster
#'
#' This function computes the log2 fold change (log2FC) of gene expression
#' for each gene between samples within a cluster and all other samples.
#'
#' It returns a matrix where each column corresponds to a cluster and each row
#' to a gene, with values representing the log2 fold change of expression in that
#' cluster versus the rest.
#'
#' @param expr A numeric matrix or data frame of normalized expression data
#'             (genes × samples). Rows are genes, columns are samples.
#' @param clusters A named or unnamed vector indicating the cluster assignment
#'                 of each sample (length equal to number of columns in `expr`).
#'
#' @return A matrix of log2 fold changes (genes × clusters).
#'         Each value is computed as log2(mean in cluster / mean in other clusters),
#'         with a small pseudocount added to avoid log of zero.
#'
#' @examples
#' \dontrun{
#' kmeans <- clustering(expr, k = 10, method = "kmeans")
#' log2fc <- compute_log2fc(expr, kmeans$clusters)
#' }
#'
#' @export
compute_log2fc <- function(expr, clusters) {
  if (!is.matrix(expr) && !is.data.frame(expr)) {
    stop("'expr' must be a numeric matrix or data frame.")
  }
  expr <- as.matrix(expr)
  if (!is.numeric(expr)) {
    stop("'expr' must contain numeric values.")
  }
  if (length(clusters) != ncol(expr)) {
    stop("Length of 'clusters' must equal the number of samples (columns) in 'expr'.")
  }

  unique_clusters <- sort(unique(clusters))
  log2fc_list <- list()

  for (cluster in unique_clusters) {
    in_cluster <- clusters == cluster
    out_cluster <- clusters != cluster

    mean_in <- rowMeans(expr[, in_cluster, drop = FALSE])
    mean_out <- rowMeans(expr[, out_cluster, drop = FALSE])

    # Add pseudocount to avoid division by zero or log(0)
    log2fc <- log2((mean_in + 1e-6) / (mean_out + 1e-6))
    log2fc_list[[as.character(cluster)]] <- log2fc
  }

  log2fc_matrix <- do.call(cbind, log2fc_list)
  colnames(log2fc_matrix) <- paste0("Cluster_", unique_clusters)
  rownames(log2fc_matrix) <- rownames(expr)

  return(log2fc_matrix)
}
