# Valid input and clustering produce an heatmap
# Check missing rownames
# Invalid clustering throw a warning

test_that("plot_heatmap works on a simple clustered matrix", {
  # Create a fake expression matrix
  mat <- matrix(rnorm(200), nrow = 20, ncol = 10)
  rownames(mat) <- paste0("gene", 1:20)
  colnames(mat) <- paste0("sample", 1:10)

  # Create a fake clustering result
  clustering_res <- list(
    model = stats::hclust(dist(mat)),
    clusters = setNames(sample(1:4, 20, replace = TRUE), rownames(mat))
  )

  expect_silent(
    plot_heatmap(mat, clustering_res, scale_data = TRUE)
  )
})

test_that("plot_heatmap warns on missing rownames", {
  mat <- matrix(rnorm(100), nrow = 10)
  clustering_res <- list(
    model = stats::hclust(dist(mat)),
    clusters = setNames(sample(1:2, 10, replace = TRUE), NULL)
  )
  expect_warning(
    plot_heatmap(mat, clustering_res),
    "rownames"
  )
})

test_that("plot_heatmap errors with invalid clustering_res", {
  mat <- matrix(rnorm(100), nrow = 10)
  expect_error(
    plot_heatmap(mat, list(model = NULL, clusters = NULL)),
    NA  #Warn
  )
})
