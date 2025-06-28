# Valid input and output
# Unvalid input throe an error

test_that("plot_gsea returns a ggplot object for a valid gseaResult", {
  # Create a fake synthetic gene ranking
  gene_ranking <- rnorm(200)
  names(gene_ranking) <- head(keys(org.Hs.eg.db, keytype = "SYMBOL"), 200)

  # Convert names to ENTREZID
  entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = names(gene_ranking),
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  entrez <- na.omit(entrez)

  # Replace names with ENTREZ to match clusterProfiler
  gene_ranking <- gene_ranking[names(entrez)]
  names(gene_ranking) <- entrez

  # Run a simple gsea
  gsea_res <- clusterProfiler::gseGO(
    geneList = sort(gene_ranking, decreasing = TRUE),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 300,
    pvalueCutoff = 0.5
  )

  # Test plot_gsea
  p <- plot_gsea(gsea_res, top_n = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plot_gsea errors on invalid input", {
  expect_error(plot_gsea(list()), "data.frame")
})
