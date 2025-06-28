#' Perform GO Enrichment Analysis for Gene Clusters
#'
#' This function performs Gene Ontology (GO) enrichment analysis on sets of genes grouped by clusters.
#' It uses the clusterProfiler package to identify enriched GO terms within each cluster. The function
#' returns one list for each cluster.
#'
#' @param gene_clusters Named vector or factor indicating cluster assignment for each gene.
#'                      The names should be gene identifiers (e.g., gene symbols or Ensembl IDs).
#' @param ontology Character string specifying which GO ontology to use.
#'                 Options are "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component).
#'                 Default is "BP".
#' @param organism_db AnnotationDbi organism database object to map gene IDs and perform enrichment.
#'                   Default is org.Hs.eg.db (human).
#' @param keyType Character string specifying the type of gene IDs used in gene_clusters.
#'                Common options include "SYMBOL" and "ENSEMBL". Default is "SYMBOL".
#'
#' @return A named list of enrichResult objects, one per cluster.
#'         Clusters with fewer than 10 mapped genes are skipped with a message.
#'
#' @importFrom clusterProfiler enrichGO
#' @importFrom AnnotationDbi mapIds
#' @export
#'
#' @examples
#' \dontrun{
#' # Suppose kmeans$clusters is a named vector with gene cluster assignments
#' go_kemeans <- go_enrichment(kmeans$clusters)
#' # Access GO enrichment for cluster "1"
#' go_kmeans[["1"]]
#' }

go_enrichment <- function(gene_clusters, organism_db) {
  if (is.null(names(gene_clusters)) || any(names(gene_clusters) == "")) {
    stop("gene_clusters must be a named vector with gene names.")
  }

  valid_genes <- keys(organism_db, keytype = "SYMBOL")
  genes_filtered <- intersect(names(gene_clusters), valid_genes)

  if (length(genes_filtered) == 0) {
    stop("No valid gene symbols found in gene_clusters for enrichment.")
  }

  gene_clusters_filtered <- gene_clusters[genes_filtered]
  enrichment_results <- list()

  for (clust in unique(gene_clusters_filtered)) {
    genes_in_cluster <- names(gene_clusters_filtered)[gene_clusters_filtered == clust]

    if (length(genes_in_cluster) < 10) {
      message(paste("Cluster", clust, "skipped due to insufficient genes"))
      next
    }

    entrez_ids <- AnnotationDbi::mapIds(
      organism_db,
      keys = genes_in_cluster,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
    entrez_ids <- na.omit(entrez_ids)

    if (length(entrez_ids) == 0) {
      warning(paste("No Entrez IDs found for cluster", clust, "- skipping enrichment."))
      next
    }

    res <- clusterProfiler::enrichGO(
      gene = entrez_ids,
      OrgDb = organism_db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )

    if (!is.null(res) && inherits(res, "enrichResult")) {
      enrichment_results[[as.character(clust)]] <- res
    } else {
      message(paste("No enrichment results for cluster", clust))
    }
  }

  return(enrichment_results)
}
