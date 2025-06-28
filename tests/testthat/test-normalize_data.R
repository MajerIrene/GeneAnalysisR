# Check all method works correctly
# Check invalid method throw error

test_that("normalize_data log2 works correctly", {
  mat <- matrix(0:5, nrow = 2)
  rownames(mat) <- c("gene1", "gene2")
  colnames(mat) <- c("s1", "s2", "s3")

  norm <- normalize_data(mat, method = "log2")

  expect_equal(norm["gene1","s1"], log2(0 + 1))
  expect_equal(norm["gene2","s3"], log2(5 + 1))
})

test_that("normalize_data zscore works correctly", {
  mat <- matrix(c(1,2,3,4,5,6), nrow = 2)
  rownames(mat) <- c("gene1", "gene2")
  colnames(mat) <- c("s1", "s2", "s3")

  norm <- normalize_data(mat, method = "zscore")

  # z-score of each gene across its own row, maintain gene names because test fail otherwise
  expect_equal(rowMeans(norm), setNames(c(0,0), c("gene1","gene2")), tolerance = 1e-8)
})

test_that("normalize_data quantile normalization preserves shapes", {
  mat <- matrix(rnorm(20), nrow = 4)
  rownames(mat) <- paste0("g", 1:4)
  colnames(mat) <- paste0("s", 1:5)

  norm <- normalize_data(mat, method = "quantile")

  expect_equal(dim(norm), dim(mat))
  expect_equal(rownames(norm), rownames(mat))
  expect_equal(colnames(norm), colnames(mat))
})

test_that("normalize_data errors on invalid method", {
  mat <- matrix(rnorm(10), nrow = 2)
  expect_error(normalize_data(mat, method = "invalid"), "should be one of")
})
