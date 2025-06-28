#' Plot GSEA Results Dotplot for Two Clusters Using a Combined Ranking
#'
#' This function visualizes Gene Set Enrichment Analysis (GSEA) results obtained
#' from a combined gene ranking of two clusters. It shows enriched pathways
#' with positive NES (cluster 1) and negative NES (cluster 2) in a single dotplot.
#'
#' Top pathways are selected by NES, regardless of statistical significance.
#' Significant terms (adjusted p < 0.05) are highlighted with stronger opacity.
#'
#' @param gsea_result An object of class gseaResult returned by gsea_enrichment.
#'
#' @param top_n Integer indicating how many top pathways to display per cluster
#'              based on NES magnitude. Default is 10.
#'
#' @return A plot displaying the dotplot of enriched pathways for both clusters.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' plot_gsea(gsea_kmeans, top_n = 10)
#' }

plot_gsea <- function(gsea_result, top_n = 10) {
  # Check added after test
  df <- tryCatch(as.data.frame(gsea_result), error = function(e) NULL)

  if (is.null(df) || !"NES" %in% colnames(df)) {
    stop("Input must be convertible to a data.frame containing an 'NES' column")
  }

  # Separate by NES direction (cluster1 positive, cluster2 negative)
  pos <- df[df$NES > 0, ]
  neg <- df[df$NES < 0, ]

  # Select top N by NES for both clusters
  top_pos <- head(pos[order(-pos$NES), ], n = top_n)
  top_neg <- head(neg[order(neg$NES), ], n = top_n)

  # Combine top positive and top negative
  combined <- rbind(top_pos, top_neg)
  combined$Cluster <- ifelse(combined$NES > 0, "Cluster 1 (positive NES)", "Cluster 2 (negative NES)")
  combined$Significant <- combined$p.adjust < 0.05
  combined$Description <- factor(combined$Description, levels = rev(unique(combined$Description)))

  # Plot
  ggplot(combined, aes(x = NES, y = Description, color = Cluster, size = -log10(p.adjust), alpha = Significant)) +
    geom_point() +
    scale_color_manual(values = c("Cluster 1 (positive NES)" = "#1b9e77", "Cluster 2 (negative NES)" = "#d95f02")) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
    labs(
      title = "GSEA Dotplot: Top Enriched Pathways in Two Clusters",
      x = "Normalized Enrichment Score (NES)",
      y = "Pathway",
      size = "-log10(adjusted p-value)",
      color = "Cluster",
      alpha = "Significant (p.adj < 0.05)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}
