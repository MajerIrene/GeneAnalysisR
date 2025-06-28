# Valid input and output
# Length mismatch in clusters

test_that("plot_pca works with valid input", {
  # simulate expression data
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(mat) <- paste0("sample", 1:10)
  colnames(mat) <- paste0("gene", 1:10)

  # simulate clustering vector
  clusters <- sample(1:2, 10, replace = TRUE)

  # check that the function returns a ggplot object without error
  p <- plot_pca(mat, clusters)
  expect_s3_class(p, "ggplot")
})

test_that("plot_pca errors with wrong cluster length", {
  mat <- matrix(rnorm(100), nrow = 10)
  wrong_clusters <- sample(1:2, 5, replace = TRUE) # length mismatch
  expect_error(plot_pca(mat, wrong_clusters))
})
