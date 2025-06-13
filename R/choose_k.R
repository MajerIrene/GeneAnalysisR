#' Determine the optimal number of clusters using elbow or silhouette method
#'
#' This function helps to determine the optimal number of clusters (\code{k}) for clustering algorithms,
#' using either the elbow method (based on within-cluster sum of squares) or the silhouette method
#' (based on average silhouette width).
#'
#' @param data A numeric data frame or matrix, with genes as rows and samples as columns.
#' @param method Character string specifying the method to use. Options are \code{"elbow"} (default)
#'               or \code{"silhouette"}.
#' @param max_k Integer indicating the maximum number of clusters to evaluate. Default is 10.
#' @param scale_data Logical, whether to scale the data by gene (row-wise z-score) before clustering. Default is \code{TRUE}.
#'
#' @return A plot is produced showing either:
#' \itemize{
#'   \item Total within-cluster sum of squares (elbow method)
#'   \item Average silhouette width (silhouette method)
#' }
#' The function is used for visual inspection and does not return a value.
#'
#' @importFrom stats kmeans dist
#' @importFrom cluster silhouette
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with gene expression matrix
#' choose_k(normalized, method = "elbow", max_k = 10)
#' choose_k(normalized, method = "silhouette", max_k = 10)
#' }
#'
choose_k <- function(data, method = "elbow", max_k = 10, scale_data = TRUE) {
  data_mat <- as.matrix(data)
  if (scale_data) {
    data_mat <- t(scale(t(data_mat)))  # scale by gene
  }

  if (method == "elbow") {
    wss <- numeric(max_k)
    for (k in 1:max_k) {
      set.seed(123)
      wss[k] <- kmeans(data_mat, centers = k)$tot.withinss
    }
    plot(1:max_k, wss, type = "b", pch = 19, col = "skyblue",
         xlab = "Number of clusters (k)",
         ylab = "Total within-cluster sum of squares",
         main = "Elbow Method")
  }

  else if (method == "silhouette") {
    if (max_k < 2) stop("Silhouette method requires max_k >= 2.")
    avg_sil <- numeric(max_k)
    for (k in 2:max_k) {
      km <- kmeans(data_mat, centers = k)
      sil <- silhouette(km$cluster, dist(data_mat))
      avg_sil[k] <- mean(sil[, 3])
    }
    plot(2:max_k, avg_sil[2:max_k], type = "b", pch = 19, col = "skyblue",
         xlab = "Number of clusters (k)",
         ylab = "Average silhouette width",
         main = "Silhouette Method")
  }

  else {
    stop("Unsupported method. Choose either 'elbow' or 'silhouette'.")
  }
}
