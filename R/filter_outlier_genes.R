#' Filter Outlier Genes Based on Expression and Variance Metrics
#'
#' This function performs filtering of genes in a gene expression matrix to remove
#' potential outlier. It filters genes based on minimum mean expression, minimum
#' variance, and extreme Z-score values of mean and variance. The function return
#' a filtered expression matrix without outlier.
#' Optional diagnostic plots visualize the distributions before filtering.
#'
#' @param expr_matrix A numeric matrix or data frame of gene expression values. Rows
#'             correspond to genes and columns to samples.
#'             It is expected that data is normalized (e.g., log-transformed).
#' @param min_mean_expr Numeric scalar. Minimum mean expression threshold.
#'                      Genes with mean expression below this are removed. Default is 1.
#' @param min_var Numeric scalar. Minimum variance threshold. Genes with variance
#'                below this are removed. Default is 0.5.
#' @param z_cutoff Numeric scalar. Z-score cutoff for detecting outliers in mean
#'                expression and variance. Genes with Z-scores outside z_cutoff
#'                are removed. Default is 3.
#' @param plot_report Logical. If TRUE, diagnostic histograms are plotted. Default is TRUE.
#'
#' @return A filtered gene expression matrix containing only genes passing all filtering criteria.
#'
#' @examples
#' \dontrun{
#' filtered_data <- filter_outlier_genes(normalized, min_mean_expr = 1, min_var = 0.5, z_cutoff = 3)
#' }
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @export
filter_outlier_genes <- function(data,
                                 min_mean_expr = 1,
                                 min_var = 0.5,
                                 z_cutoff = 3,
                                 plot_report = TRUE) {

  # Controllo corretto per rownames
  if (is.null(rownames(data))) {
    stop("Input expression matrix must have row names (gene names).")
  }

  z_score <- function(x) (x - mean(x)) / sd(x)

  gene_means <- rowMeans(data)
  gene_vars <- apply(data, 1, var)

  z_means <- z_score(gene_means)
  z_vars <- z_score(gene_vars)

  keep_expr <- gene_means > min_mean_expr
  keep_var <- gene_vars > min_var
  keep_z <- (abs(z_means) < z_cutoff) & (abs(z_vars) < z_cutoff)

  keep <- keep_expr & keep_var & keep_z

  data_filtered <- data[keep, , drop = FALSE]

  if (plot_report) {
    p1 <- ggplot(data.frame(value = gene_means), aes(x = value)) +
      geom_histogram(bins = 50, fill = "skyblue", color = "black") +
      geom_vline(xintercept = min_mean_expr, color = "red", linetype = "dashed") +
      ggtitle("Distribution of Mean Gene Expression")

    p2 <- ggplot(data.frame(value = gene_vars), aes(x = value)) +
      geom_histogram(bins = 50, fill = "lightblue", color = "black") +
      geom_vline(xintercept = min_var, color = "red", linetype = "dashed") +
      ggtitle("Distribution of Gene Variance")

    p3 <- ggplot(data.frame(z_means = z_means), aes(x = z_means)) +
      geom_histogram(bins = 50, fill = "deepskyblue", color = "black") +
      geom_vline(xintercept = c(-z_cutoff, z_cutoff), color = "red", linetype = "dashed") +
      ggtitle("Z-score of Mean Expression")

    p4 <- ggplot(data.frame(z_vars = z_vars), aes(x = z_vars)) +
      geom_histogram(bins = 50, fill = "dodgerblue", color = "black") +
      geom_vline(xintercept = c(-z_cutoff, z_cutoff), color = "red", linetype = "dashed") +
      ggtitle("Z-score of Variance")

    gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

    message(sprintf("Total genes: %d", nrow(data)))
    message(sprintf("Filtered genes: %d (%.2f%% retained)", nrow(data_filtered), 100 * nrow(data_filtered) / nrow(data)))
  }

  return(data_filtered)
}
