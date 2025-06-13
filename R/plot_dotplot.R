#' Plot a Dotplot for a Single GO Enrichment Result
#'
#' This function creates a dotplot from a single GO enrichment result (from `enrichGO`, `enrichKEGG`, etc.).
#'
#' @param enrichment_res An object of class `enrichResult`, typically the result of `enrichGO()`.
#' @param title A character string specifying the plot title.
#'
#' @return A ggplot object representing the dotplot of enrichment results.
#'
#' @importFrom enrichplot dotplot
#' @importFrom ggplot2 ggtitle
#' @export
#'
#' @examples
#' \dontrun{
#' plot_dotplot(go_kmeans[[1]], title = "GO Enrichment - KMeans cluster 1")
#' }
plot_dotplot <- function(enrichment_res, title = "GO Enrichment Dotplot") {
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("The 'enrichplot' package is required. Please install it using BiocManager::install('enrichplot').")
  }
  # Can't pass any argument function because of errors
  p <- enrichplot::dotplot(enrichment_res) + ggplot2::ggtitle(title)
  return(p)
}
