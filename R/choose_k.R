#' Determine the optimal number of clusters using elbow and silhouette method
#'
#' This function helps to determine the optimal number of clusters (\code{k}) for clustering algorithms,
#' using both the elbow method or the silhouette method.
#'
#' @param expr_data Numeric matrix of gene expression (genes x samples).
#' @param method Character string specifying the method to use. Options are \code{"elbow"} (default)
#'               or \code{"silhouette"}.
#' @param max_k Integer indicating the maximum number of clusters to evaluate. Default is 10.
#'
#' @return A plot is produced showing either:
#' \itemize{
#'   \item Total within-cluster sum of squares (elbow method)
#'   \item Average silhouette width (silhouette method)
#' }
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

choose_k <- function(expr_data, method = "elbow", max_k = 10, scale_data = TRUE) {
  data_mat <- as.matrix(expr_data)

  # Elbow plot
  if (method == "elbow") {
    wss <- numeric(max_k)
    for (k in 1:max_k) {
      set.seed(123)
      # Using kmeans to estimate elbow plot
      wss[k] <- kmeans(data_mat, centers = k)$tot.withinss
    }
    plot(1:max_k, wss, type = "b", pch = 19, col = "skyblue",
         xlab = "Number of clusters (k)",
         ylab = "Total within-cluster sum of squares",
         main = "Elbow Method")
  }

  # Silhouette
  else if (method == "silhouette") {
    if (max_k < 2) stop("Silhouette method requires max_k >= 2.")
    avg_sil <- numeric(max_k)
    for (k in 2:max_k) {
      # Using kmeans to estimate the silhouette
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
