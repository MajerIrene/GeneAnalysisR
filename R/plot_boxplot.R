#' Plot Boxplot of Expression Data
#'
#' This function plots a boxplot of gene expression values across samples,
#' either for raw or normalized data.
#'
#' @param expr_data A numeric matrix or data frame of expression data (genes x samples).
#' @param title A character string for the plot title (e.g., "Raw Data" or "Normalized Data").
#'
#' @return A ggplot2 object with the boxplot.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw theme element_text labs
#' @importFrom reshape2 melt
#' @export
#'
#' @examples
#' \dontrun{
#'   plot_boxplot(expr_data = raw, title = "Raw Expression")
#'   plot_boxplot(expr_data = normalized, title = "Normalized Expression")
#' }
plot_boxplot <- function(expr_data, title = "Expression Data") {
  long_data <- melt(expr_data)
  colnames(long_data) <- c("Gene", "Sample", "Expression")

  ggplot(long_data, aes(x = Sample, y = Expression)) +
    geom_boxplot(outlier.size = 0.5, fill = "skyblue") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, y = "Expression Level", x = "Sample")
}
