# Check return type of the enrichment
# Cluster with less than 10 genes are skipped
# Cluster with more than 10 genes are enriched
# Cluster without names throw an error

test_that("go_enrichment returns a list", {
  gene_clusters <- c(
    rep("1", 12),
    rep("2", 8)
  )
  real_genes <- c(
    "TP53", "BRCA1", "EGFR", "MYC", "CDK2", "MDM2", "PTEN", "AKT1",
    "CCND1", "CDKN1A", "BRAF", "KRAS",
    "MTOR", "MAPK1", "PIK3CA", "RB1", "CCNE1", "SMAD4", "FGFR1", "JAK2"
  )
  names(gene_clusters) <- real_genes

  res <- go_enrichment(gene_clusters, organism_db = org.Hs.eg.db)
  expect_type(res, "list")
})

test_that("clusters with fewer than 10 genes are skipped with message", {
  gene_clusters <- c(
    rep("1", 5),
    rep("2", 3)
  )
  real_genes <- c(
    "TP53", "BRCA1", "EGFR", "MYC", "CDK2",
    "MDM2", "PTEN", "AKT1"
  )
  names(gene_clusters) <- real_genes

  expect_message(
    go_enrichment(gene_clusters, organism_db = org.Hs.eg.db),
    "skipped due to insufficient genes"
  )
})

test_that("result contains enrichResult objects for clusters with >=10 genes", {
  gene_clusters <- c(
    rep("1", 15),
    rep("2", 5)
  )
  real_genes <- c(
    "TP53", "BRCA1", "EGFR", "MYC", "CDK2",
    "MDM2", "PTEN", "AKT1", "CCND1", "CDKN1A",
    "BRAF", "KRAS", "RB1", "CCNE1", "SMAD4",   # 15 genes for cluster 1
    "MAPK1", "PIK3CA", "FGFR1", "JAK2", "MTOR"  # 5 genes for cluster 2
  )
  names(gene_clusters) <- real_genes

  res <- go_enrichment(gene_clusters, organism_db = org.Hs.eg.db)
  expect_true("1" %in% names(res))
  expect_true(inherits(res[["1"]], "enrichResult"))
  expect_false("2" %in% names(res))  # cluster 2 too small, should be skipped
})

test_that("error if gene_clusters has no names", {
  gene_clusters <- factor(rep("1", 12))
  expect_error(
    go_enrichment(gene_clusters, organism_db = org.Hs.eg.db),
    "must be a named vector"
  )
})
