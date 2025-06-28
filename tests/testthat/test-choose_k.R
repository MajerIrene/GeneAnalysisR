# Elbow no error
# Silhouette no error
# Invalid method
# k<2 silhouette
# TBD

test_that("choose_k elbow method runs without errors", {
  mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  rownames(mat) <- paste0("gene", 1:100)

  expect_silent(choose_k(mat, method = "elbow", max_k = 5))
})

test_that("choose_k silhouette method runs without errors", {
  mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  rownames(mat) <- paste0("gene", 1:100)

  expect_silent(choose_k(mat, method = "silhouette", max_k = 5))
})

test_that("choose_k stops with unsupported method", {
  mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  rownames(mat) <- paste0("gene", 1:100)

  expect_error(choose_k(mat, method = "invalid"), "Unsupported method")
})

test_that("choose_k stops with too small max_k for silhouette", {
  mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  rownames(mat) <- paste0("gene", 1:100)

  expect_error(choose_k(mat, method = "silhouette", max_k = 1), "requires max_k >= 2")
})
