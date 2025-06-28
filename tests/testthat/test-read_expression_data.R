# CHeck CSV and TSV
# File is missing
# Unsupported file format
# Non numeric data
# No rownames

test_that("read_expression_data reads CSV correctly", {
  # Crea un file temporaneo CSV
  tmp <- tempfile(fileext = ".csv")
  write.csv(data.frame(Gene1 = 1:3, Gene2 = 4:6), tmp, row.names = TRUE)

  mat <- read_expression_data(tmp)

  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(3, 2))
  expect_equal(rownames(mat), as.character(1:3))
  expect_equal(colnames(mat), c("Gene1", "Gene2"))
  expect_equal(mat[1,1], 1)
})

test_that("read_expression_data reads TSV correctly", {
  tmp <- tempfile(fileext = ".tsv")
  write.table(data.frame(Gene1 = 1:3, Gene2 = 4:6), tmp, sep = "\t", row.names = TRUE, quote = FALSE)

  mat <- read_expression_data(tmp)

  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(3, 2))
})

test_that("read_expression_data errors if file missing", {
  expect_error(read_expression_data("nonexistentfile.csv"), "File does not exist")
})

test_that("read_expression_data errors on unsupported extension", {
  tmp <- tempfile(fileext = ".xlsx")
  file.create(tmp)
  expect_error(read_expression_data(tmp), "Unsupported file extension")
})

test_that("read_expression_data errors on non-numeric data", {
  tmp <- tempfile(fileext = ".csv")
  df <- data.frame(Gene1 = c("a", "b", "c"), Gene2 = c(1, 2, 3))
  write.csv(df, tmp, row.names = TRUE)

  expect_error(read_expression_data(tmp), "Non-numeric values detected")
})

test_that("read_expression_data handles no row names", {
  tmp <- tempfile(fileext = ".csv")
  write.csv(data.frame(Gene1 = 1:3, Gene2 = 4:6), tmp, row.names = FALSE)

  mat <- read_expression_data(tmp, row_names = FALSE)

  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 3)
})
