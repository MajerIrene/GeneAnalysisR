# Correct input produces correct output (Check better)
# Model is missing
# Model is not hierarchical
# TBD


test_that("plot_dendrogram produces a plot for valid clustering result", {
  # Simulated hierarchical clustering
  mat <- matrix(rnorm(40), nrow = 10)
  rownames(mat) <- paste0("gene", 1:10)
  dist_mat <- dist(mat)
  hc <- hclust(dist_mat)

  clusters <- cutree(hc, k = 3)

  res <- list(
    model = hc,
    clusters = clusters
  )

  expect_silent(plot_dendrogram(res))
})


test_that("plot_dendrogram errors if list lacks model or clusters", {
  bad1 <- list(model = hclust(dist(matrix(rnorm(10), nrow = 5))))
  bad2 <- list(clusters = rep(1, 5))

  expect_error(plot_dendrogram(bad1), "clusters")
  expect_error(plot_dendrogram(bad2), "model")
})

test_that("plot_dendrogram errors if model is not of class hclust", {
  fake_model <- list(dummy = TRUE)
  res <- list(model = fake_model, clusters = rep(1, 5))

  expect_error(plot_dendrogram(res), "Model is not hierarchical")
})
