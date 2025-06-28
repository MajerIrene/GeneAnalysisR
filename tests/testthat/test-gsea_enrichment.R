# Output check
# Empty input
# Check different ontology
# Throw error on missing names
# TBD

test_that("gsea_enrichment returns gseaResult", {
  # Use real gene symbols from OrgDb keys
  real_genes <- sample(keys(org.Hs.eg.db, keytype = "SYMBOL"), 100)

  gene_ranking <- rnorm(100)
  names(gene_ranking) <- real_genes

  res <- gsea_enrichment(
    gene_ranking = gene_ranking,
    OrgDb = org.Hs.eg.db
  )

  expect_s4_class(res, "gseaResult")
})

test_that("gsea_enrichment errors on empty ranking", {
  gene_ranking <- numeric(0)
  names(gene_ranking) <- character(0)

  expect_error(
    gsea_enrichment(
      gene_ranking = gene_ranking,
      OrgDb = org.Hs.eg.db
    ),
    "gene_ranking must be a non-empty named numeric vector"
  )
})

test_that("gsea_enrichment works with different ontologies", {
  real_genes <- sample(keys(org.Hs.eg.db, keytype = "SYMBOL"), 100)

  gene_ranking <- rnorm(100)
  names(gene_ranking) <- real_genes

  res <- gsea_enrichment(
    gene_ranking = gene_ranking,
    OrgDb = org.Hs.eg.db,
    ont = "MF"
  )

  expect_s4_class(res, "gseaResult")
})

test_that("gsea_enrichment errors on missing names", {
  gene_ranking <- rnorm(50)
  names(gene_ranking) <- NULL

  expect_error(
    gsea_enrichment(
      gene_ranking = gene_ranking,
      OrgDb = org.Hs.eg.db
    ),
    "gene_ranking must be a non-empty named numeric vector"
  )
})
