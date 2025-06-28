#' Compute gene ranking for GSEA based on correlation to cluster centroids
#'
#' This function calculates the correlation of each gene's expression profile
#' within the centroid profiles of two gene clusters, then computes the difference
#' in correlation to rank genes for GSEA. Positive numbers are for cluster1 and
#' negative for cluster2
#'
#' @param expr_mat Numeric matrix of gene expression (genes x samples).
#' @param gene_clusters Named vector of cluster assignments for genes. Results of
#' clustering function
#' @param cluster1 Cluster number or name for the first cluster.
#' @param cluster2 Cluster number or name for the second cluster.
#'
#' @return Named numeric vector of gene ranking scores (correlation difference),
#' sorted in decreasing order.
#'
#' @export
#' @examples
#' \dontrun{
#' ranking <- compute_gene_ranking(normalized, kmeans, 1, 2)
#' }
compute_gene_ranking <- function(expr_mat, gene_clusters, cluster1, cluster2) {
  # Adding check for test
  if (is.null(rownames(expr_mat))) stop("expr_mat must have row names (gene names).")
  if (is.null(names(gene_clusters))) stop("gene_clusters must be a named vector with gene names.")

  genes <- rownames(expr_mat)
  if (!all(names(gene_clusters) %in% genes)) stop("All gene_clusters names must be in expr_mat rownames.")

  genes_c1 <- names(gene_clusters)[gene_clusters == cluster1]
  genes_c2 <- names(gene_clusters)[gene_clusters == cluster2]

  # Compute centroid profiles if cluster has genes, else NA vector -> Test failed before so better to add this
  centroid1 <- if (length(genes_c1) > 0) {
    colMeans(expr_mat[genes_c1, , drop = FALSE])
  } else {
    rep(NA_real_, ncol(expr_mat))
  }
  centroid2 <- if (length(genes_c2) > 0) {
    colMeans(expr_mat[genes_c2, , drop = FALSE])
  } else {
    rep(NA_real_, ncol(expr_mat))
  }

  # Safe corr to deal with NA
  safe_cor <- function(gene_expr, centroid) {
    if (all(is.na(centroid))) return(NA_real_)
    cor(gene_expr, centroid, use = "complete.obs")
  }

  # For each gene compute correlation to each centroid
  cor_to_c1 <- apply(expr_mat, 1, function(gene_expr) safe_cor(gene_expr, centroid1))
  cor_to_c2 <- apply(expr_mat, 1, function(gene_expr) safe_cor(gene_expr, centroid2))

  # If one centroid is NA for all, treat correlation as 0 for that cluster
  cor_to_c1[is.na(cor_to_c1)] <- 0
  cor_to_c2[is.na(cor_to_c2)] <- 0

  ranking_score <- cor_to_c1 - cor_to_c2

  # Return ranking scores named by genes
  ranking_score <- ranking_score[!is.na(ranking_score)]
  ranking_score <- sort(ranking_score, decreasing = TRUE) #ordered decreasing
  return(ranking_score)
}
