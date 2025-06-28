#' Plot a Dotplot for a Single GO Enrichment Result
#'
#' This function creates a dotplot from a single GO enrichment result.
#'
#' @param enrichment_res An object of class `enrichResult`, that is the result of `enrichGO()`.
#' @param showCategory Number of enriched term to show in the plot. Default is 10.
#' @param title A character string specifying the plot title.
#'
#' @return Plots the dotplot of enrichment results.
#'
#' @importFrom enrichplot dotplot
#' @importFrom ggplot2 ggtitle
#' @export
#'
#' @examples
#' \dontrun{
#' plot_dotplot(go_kmeans[[1]], title = "GO Enrichment - KMeans cluster 1", showCategory = 15)
#' }
plot_dotplot <- function(enrichment_res, title, showCategory = 10) {
  # Adding check after test
  valid_classes <- c("enrichResult", "gseaResult", "compareClusterResult")

  if (!inherits(enrichment_res, valid_classes)) {
    stop("Input must be an enrichment result object of class enrichResult, gseaResult, or compareClusterResult")
  }

  p <- dotplot(enrichment_res, showCategory = showCategory) + ggtitle(title)
  return(p)
}
