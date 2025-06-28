# Clustering works with different methods
# K=0 invalid
# k>sample invalid
# Non numeric data
# Missing rownames
# TBD


test_that("clustering works with kmeans", {
  set.seed(42)
  expr_data <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5)
  rownames(expr_data) <- paste0("gene", 1:50)

  res <- clustering(expr_data, k = 3, method = "kmeans")

  expect_type(res, "list")
  expect_named(res, c("clusters", "model"))
  expect_length(res$clusters, 50)
  expect_true(all(names(res$clusters) == rownames(expr_data)))
  expect_s3_class(res$model, "kmeans")
})

test_that("clustering works with hierarchical", {
  set.seed(42)
  expr_data <- matrix(rnorm(30 * 5), nrow = 30, ncol = 5)
  rownames(expr_data) <- paste0("gene", 1:30)

  res <- clustering(expr_data, k = 4, method = "hierarchical")

  expect_type(res, "list")
  expect_named(res, c("clusters", "model"))
  expect_length(res$clusters, 30)
  expect_s3_class(res$model, "hclust")
})

test_that("clustering works with pam", {
  skip_if_not_installed("cluster")

  set.seed(42)
  expr_data <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  rownames(expr_data) <- paste0("gene", 1:20)

  res <- clustering(expr_data, k = 3, method = "pam")

  expect_type(res, "list")
  expect_named(res, c("clusters", "model"))
  expect_length(res$clusters, 20)
  expect_s3_class(res$model, "pam")
})

test_that("clustering fails with invalid k", {
  expr_data <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  rownames(expr_data) <- paste0("gene", 1:10)

  expect_error(clustering(expr_data, k = 0, method = "kmeans"), "must be >1")
  expect_error(clustering(expr_data, k = 20, method = "kmeans"), "must be >1")
})

test_that("clustering fails with missing rownames", {
  expr_data <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  # no rownames
  expect_error(clustering(expr_data, k = 3, method = "kmeans"), "must have row names")
})

test_that("clustering fails with non-numeric data", {
  expr_data <- matrix(letters[1:20], nrow = 10, ncol = 2)
  rownames(expr_data) <- paste0("gene", 1:10)
  expect_error(clustering(expr_data, k = 2, method = "kmeans"), "must be numeric")
})


