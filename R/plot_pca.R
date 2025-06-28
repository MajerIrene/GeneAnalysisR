#' Plot PCA colored by clustering results
#'
#' This function visualizes the first two principal components of a dataset,
#' coloring the points based on cluster assignments.
#'
#' @param expr_data A numeric matrix or data frame used for clustering (rows = samples).
#' @param clusters A vector of cluster assignments (same length as number of rows in data).
#' @return NULL. Produces a PCA plot.
#' @import stats
#' @import graphics
#' @export
#' @examples
#' \dontrun{
#' plot_pca(normalized, kmeans$clusters)
#' }
plot_pca <- function(expr_data, clusters) {
  if (length(clusters) != nrow(expr_data)) {
    stop("Length of clusters vector does not match number of samples (rows) in expr_data.")
  }

  data_numeric <- as.matrix(expr_data)

  pca <- prcomp(data_numeric, scale. = TRUE)
  pc_data <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Cluster = factor(clusters)
  )

  explained_var <- round(100 * summary(pca)$importance[2, 1:2], 1)

  colors <- c("#084594", "#bdd7e7")

  ggplot(pc_data, aes(x = PC1, y = PC2, color = Cluster)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = colors) +
    labs(
      title = "PCA Plot",
      x = paste0("PC1 (", explained_var[1], "%)"),
      y = paste0("PC2 (", explained_var[2], "%)"),
      color = "Cluster"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}
