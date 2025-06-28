# Gene under threshold are removed
# All gene pass
# Check the plot_report parameter
# Gene names missing (Makes sense?)

test_that("filter_outlier_genes filters genes based on thresholds", {
  set.seed(123)
  expr <- matrix(rnorm(50, mean = 5, sd = 2), nrow = 10)
  rownames(expr) <- paste0("Gene", 1:10)

  expr["Gene1", ] <- rep(0.1, 5)

  filtered <- filter_outlier_genes(expr, min_mean_expr = 1, min_var = 0.5, z_cutoff = 3, plot_report = FALSE)

  expect_false("Gene1" %in% rownames(filtered))
  expect_lte(nrow(filtered), nrow(expr))
})

test_that("filter_outlier_genes returns same data if all genes pass", {
  set.seed(456)
  expr <- matrix(rnorm(50, mean = 10, sd = 5), nrow = 10)
  rownames(expr) <- paste0("Gene", 1:10)

  filtered <- filter_outlier_genes(expr, min_mean_expr = 1, min_var = 0.5, z_cutoff = 3, plot_report = FALSE)

  expect_equal(nrow(filtered), nrow(expr))
  expect_equal(rownames(filtered), rownames(expr))
})

test_that("filter_outlier_genes respects plot_report argument", {
  expr <- matrix(rnorm(50, mean = 5), nrow = 10)
  rownames(expr) <- paste0("Gene", 1:10)

  expect_message(filter_outlier_genes(expr, plot_report = TRUE), "Total genes")
  expect_silent(filter_outlier_genes(expr, plot_report = FALSE))
})

test_that("filter_outlier_genes errors on missing rownames", {
  expr <- matrix(rnorm(50, mean = 5), nrow = 10)

  expect_error(filter_outlier_genes(expr), "subscript out of bounds|row names")
})
