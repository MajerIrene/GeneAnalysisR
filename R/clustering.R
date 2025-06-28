#' Perform clustering using k-means, hierarchical, or PAM method
#'
#' This function performs either k-means, hierarchical or PAM clustering over the
#' gene set (not samples). The function takes as input the expression matrix
#' (genes x samples) and returns a list containing cluster assignment and the
#' chosen model.
#'
#' @param expr_data Numeric matrix of gene expression (genes x samples).
#' @param k Number of clusters.
#' @param method Clustering method: "kmeans", "hierarchical", or "pam".
#' @param nstart Number of starts for k-means. Default is 25.
#' @param distance_method Distance metric for hierarchical or PAM. Default is euclidean.
#'
#' @return A list with:
#' \itemize{
#'  \item clusters: Named vector of cluster assignments (gene names as names).
#'  \item model: The clustering model object (kmeans, hclust, or pam).
#' }
#' @import stats
#' @importFrom cluster pam
#' @export
clustering <- function(expr_data, k, method = c("kmeans", "hierarchical", "pam"), nstart = 25, distance_method = "euclidean") {
  method <- match.arg(method)

  if (!is.matrix(expr_data) && !is.data.frame(expr_data)) stop("'expr_data' must be a matrix or data frame")
  expr_data <- as.matrix(expr_data)
  if (!is.numeric(expr_data)) stop("'expr_data' must be numeric")
  if (k <= 1 || k > nrow(expr_data)) stop("'k' must be >1 and <= number of rows")

  if (is.null(rownames(expr_data))) stop("'expr_data' must have row names for gene identifiers")

  if (method == "kmeans") {
    model <- kmeans(expr_data, centers = k, nstart = nstart)
    clusters <- model$cluster

  } else if (method == "hierarchical") {
    dist_mat <- dist(expr_data, method = distance_method)
    model <- hclust(dist_mat, method = "ward.D2")
    clusters <- cutree(model, k = k)

  } else if (method == "pam") {
    dist_mat <- dist(expr_data, method = distance_method)
    model <- pam(dist_mat, k = k)
    clusters <- model$clustering
  }

  names(clusters) <- rownames(expr_data)

  return(list(clusters = clusters, model = model))
}
