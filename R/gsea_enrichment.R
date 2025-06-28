#' Perform GSEA using clusterProfiler and customizable organism
#'
#' This function runs Gene Set Enrichment Analysis (GSEA) using the package clusterProfiler,
#' assuming the input ranking is a named numeric vector where names are gene identifiers
#' (e.g., Ensembl IDs). The user can specify the organism database, key type, ontologies,
#' minimum and maximum size of a gene set and the p-value cutoff.
#'
#' @param gene_ranking Named numeric vector of gene scores (e.g., log2FC), names are gene IDs.
#' @param OrgDb Annotation package, e.g., org.Hs.eg.db or org.Mm.eg.db. Default is org.Hs.eg.db.
#' @param keyType The type of gene IDs used in ranking (e.g., "ENSEMBL", "SYMBOL"). Default is SYMBOL.
#' @param ont Ontology for GO analysis: "BP", "MF", or "CC". Default is "BP".
#' @param minGSSize Minimum size of a gene set. Default is 10.
#' @param maxGSSize Maximum size of a gene set. Default is 500.
#' @param pvalueCutoff Adjusted p-value cutoff. Default is 0.05.
#'
#' @return A gseaResult object.
#'
#' @importFrom clusterProfiler gseGO
#' @export
#'
#' @examples
#' \dontrun{
#'   result <- gsea_enrichment(
#'     gene_ranking = ranking_kmeans,
#'     OrgDb = org.Hs.eg.db,
#'     keyType = "SYMBOL",
#'     ont = "BP"
#'   )
#' }
gsea_enrichment <- function(gene_ranking,
                                     OrgDb,
                                     keyType = "SYMBOL",
                                     ont = "BP",
                                     minGSSize = 10,
                                     maxGSSize = 500,
                                     pvalueCutoff = 0.05) {
  gene_ranking <- sort(gene_ranking, decreasing = TRUE)
  gene_ranking <- gene_ranking[!is.na(names(gene_ranking)) & names(gene_ranking) != ""]

  gsea_result <- gseGO(
    geneList = gene_ranking,
    OrgDb = OrgDb,
    keyType = keyType,
    ont = ont,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = pvalueCutoff,
    verbose = TRUE
  )

  return(gsea_result)
}
