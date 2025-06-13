#' Perform GO Enrichment Analysis for Gene Clusters
#'
#' This function performs Gene Ontology (GO) enrichment analysis on sets of genes grouped by clusters.
#' It uses the \code{clusterProfiler} package to identify enriched GO terms within each cluster.
#'
#' @param gene_clusters Named vector or factor indicating cluster assignment for each gene.
#'                      The names should be gene identifiers (e.g., gene symbols or Ensembl IDs).
#' @param ontology Character string specifying which GO ontology to use.
#'                 Options are \code{"BP"} (Biological Process), \code{"MF"} (Molecular Function), or \code{"CC"} (Cellular Component).
#'                 Default is \code{"BP"}.
#' @param organism_db AnnotationDbi organism database object to map gene IDs and perform enrichment.
#'                   Default is \code{org.Hs.eg.db::org.Hs.eg.db} (human).
#' @param keyType Character string specifying the type of gene IDs used in \code{gene_clusters}.
#'                Common options include \code{"SYMBOL"} and \code{"ENSEMBL"}. Default is \code{"SYMBOL"}.
#'
#' @return A named list of \code{enrichResult} objects (from \code{clusterProfiler}), one per cluster.
#'         Clusters with fewer than 10 mapped genes are skipped with a message.
#'
#' @importFrom clusterProfiler enrichGO
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @export
#'
#' @examples
#' \dontrun{
#' # Suppose kmeans$clusters is a named vector with gene cluster assignments
#' go_kemeans <- go_enrichment(kmeans$clusters)
#' # Access GO enrichment for cluster "1"
#' go_results[["1"]]
#' }

go_enrichment <- function(gene_clusters,
                          ontology = "BP",
                          organism_db = org.Hs.eg.db::org.Hs.eg.db,
                          keyType = "SYMBOL") {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Please install the 'clusterProfiler' package.")
  }

  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Please install the 'AnnotationDbi' package.")
  }

  results_list <- list()
  cluster_ids <- unique(gene_clusters)

  for (clust in cluster_ids) {
    gene_subset <- names(gene_clusters[gene_clusters == clust])
    gene_subset <- unique(gene_subset)

    entrez_ids <- suppressMessages(
      AnnotationDbi::mapIds(
        x = organism_db,
        keys = gene_subset,
        column = "ENTREZID",
        keytype = keyType,
        multiVals = "first"
      )
    )

    entrez_ids <- na.omit(entrez_ids)

    if (length(entrez_ids) >= 10) {
      go_result <- clusterProfiler::enrichGO(
        gene         = entrez_ids,
        OrgDb        = organism_db,
        keyType      = "ENTREZID",
        ont          = ontology,
        pAdjustMethod= "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable     = TRUE
      )
      results_list[[as.character(clust)]] <- go_result
    } else {
      message(sprintf("Cluster %s skipped: fewer than 10 valid genes.", clust))
    }
  }

  return(results_list)
}
