# Check return type as ggplot
# Boxplot uses provided title
# Check for non numeric input
# Other check for input?

test_that("plot_boxplot returns a ggplot object", {
  mat <- matrix(rnorm(20), nrow = 4)
  rownames(mat) <- paste0("gene", 1:4)
  colnames(mat) <- paste0("sample", 1:5)

  p <- plot_boxplot(mat, title = "Test Plot")

  expect_s3_class(p, "ggplot")
})

test_that("plot_boxplot uses the provided title", {
  mat <- matrix(rnorm(12), nrow = 3)
  rownames(mat) <- paste0("g", 1:3)
  colnames(mat) <- paste0("s", 1:4)

  p <- plot_boxplot(mat, title = "My Custom Title")

  expect_match(p$labels$title, "My Custom Title")
})

test_that("plot_boxplot fails on non-numeric input", {
  bad <- matrix(letters[1:12], nrow = 3)
  expect_error(plot_boxplot(bad), "numeric") # more tolerant
})

