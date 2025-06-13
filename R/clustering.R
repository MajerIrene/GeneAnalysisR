#' Perform clustering using k-means, hierarchical, or PAM method
#'
#' This function performs clustering on a numeric dataset using either k-means, hierarchical, or PAM clustering,
#' based on the user's choice. It returns the clustering results.
#'
#' @param data A numeric matrix or data frame containing the features to cluster (genes × samples).
#' @param k Integer specifying the number of clusters.
#' @param method Character string specifying the clustering method: "kmeans", "hierarchical", or "pam".
#' @param nstart Number of random starts for k-means (default 25). Ignored if method is "hierarchical" or "pam".
#' @param distance_method Distance metric to use for hierarchical or PAM clustering (default "euclidean").
#'        Options include "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski".
#'
#' @return A list containing clustering results:
#' \itemize{
#'   \item \code{clusters}: a named vector of cluster assignments.
#'   \item \code{model}: the clustering model object (kmeans, hclust, or pam).
#' }
#' @import stats
#' @importFrom cluster pam
#' @examples
#' \dontrun{
#' kmeans <- clustering(normalized, k = 10, method = "kmeans")
#' hierarchical<- clustering(normalized, k = 10, method = "hierarchical", distance_method = "manhattan")
#' pam <- clustering(normalized, k = 10, method = "pam")
#' }
#' @export
clustering <- function(data, k, method = c("kmeans", "hierarchical", "pam"),
                       nstart = 25, distance_method = "euclidean") {
  method <- match.arg(method)

  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a numeric matrix or data frame")
  }
  data <- as.matrix(data)
  if (!is.numeric(data)) {
    stop("'data' must be numeric")
  }
  if (!is.numeric(k) || length(k) != 1 || k <= 1) {
    stop("'k' must be a single integer greater than 1")
  }
  if (k > nrow(data)) {
    warning("'k' is larger than the number of rows in 'data'")
  }

  if (method == "kmeans") {
    model <- kmeans(data, centers = k, nstart = nstart)
    clusters <- model$cluster

  } else if (method == "hierarchical") {
    dist_mat <- dist(data, method = distance_method)
    model <- hclust(dist_mat, method = "ward.D2")
    clusters <- cutree(model, k = k)

  } else if (method == "pam") {
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop("Please install the 'cluster' package to use PAM.")
    }
    dist_mat <- dist(data, method = distance_method)
    model <- cluster::pam(dist_mat, k = k)
    clusters <- model$clustering
  }

  if (is.null(rownames(data))) {
    stop("Input data must have rownames corresponding to gene symbols.")
  }
  names(clusters) <- rownames(data)

  return(list(clusters = clusters, model = model))
}
