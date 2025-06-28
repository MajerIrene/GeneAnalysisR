# Output check
# Empty input
# Check different ontology
# Throw error on missing names
# TBD


test_that("gsea_enrichment returns gseaResult", {
  # create ranking for test
  gene_ranking <- rnorm(100)
  names(gene_ranking) <- paste0("GENE", seq_along(gene_ranking))

  res <- gsea_enrichment(
    gene_ranking = gene_ranking,
    OrgDb = org.Hs.eg.db
  )

  expect_s3_class(res, "gseaResult")
})

test_that("gsea_enrichment handles empty ranking", {
  gene_ranking <- numeric(0)
  names(gene_ranking) <- character(0)

  expect_error(
    gsea_enrichment(
      gene_ranking,
      OrgDb = org.Hs.eg.db
    ),
    "geneList should be a named numeric vector" # error comes from gseGO
  )
})

test_that("gsea_enrichment works with different ontologies", {
  gene_ranking <- rnorm(100)
  names(gene_ranking) <- paste0("GENE", seq_along(gene_ranking))

  res <- gsea_enrichment(
    gene_ranking,
    OrgDb = org.Hs.eg.db,
    ont = "MF"
  )

  expect_s3_class(res, "gseaResult")
})

test_that("gsea_enrichment errors on missing names", {
  gene_ranking <- rnorm(50)
  names(gene_ranking) <- NULL

  expect_error(
    gsea_enrichment(
      gene_ranking,
      OrgDb = org.Hs.eg.db
    ),
    "geneList should be a named numeric vector" # from gseGO
  )
})
