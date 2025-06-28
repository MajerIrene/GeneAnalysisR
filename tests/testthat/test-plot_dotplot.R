# Valid input and output
# Invalid input throw an error
# TBD

test_that("plot_dotplot returns a ggplot object for a valid enrichResult", {
  # Create a fake enrichResult object
  genes <- rep(1, 20)
  names(genes) <- head(keys(org.Hs.eg.db, keytype = "SYMBOL"), 20)
  entrez <- mapIds(
    org.Hs.eg.db,
    keys = names(genes),
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  entrez <- na.omit(entrez)

  enrich <- enrichGO(
    gene = entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5
  )

  p <- plot_dotplot(enrich, title = "Test Dotplot", showCategory = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plot_dotplot errors with invalid enrichment object", {
  expect_error(plot_dotplot(list(), title = "invalid"), "enrichResult")
})
