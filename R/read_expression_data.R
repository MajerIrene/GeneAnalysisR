#' Read gene expression data from file
#'
#' Reads a gene expression matrix from a CSV, TSV, or TXT file. Assumes genes are in rows and samples in columns.
#'
#' @param file_path Path to the file (CSV, TSV, or TXT).
#' @param header Logical. Whether the file contains a header row. Default is TRUE.
#' @param row_names Logical. Whether the first column contains gene names. Default is TRUE.
#'
#' @return A numeric matrix of gene expression values.
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "airway_top500.csv", package = "GeneClustR")
#' expr <- read_expression_data(file_path)

read_expression_data <- function(file_path, header = TRUE, row_names = TRUE) {
  # Check file existence
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }

  # Use different separator for different file format
  ext <- tools::file_ext(file_path)
  sep <- switch(
    tolower(ext),
    "csv" = ",",
    "tsv" = "\t",
    "txt" = "\t",
    stop("Unsupported file extension: file format must be .csv, .tsv, or .txt")
  )

  # Read file
  df <- tryCatch(
    {
      read.delim(file_path, sep = sep, header = header, row.names = if (row_names) 1 else NULL, check.names = FALSE)
    },
    error = function(e) {
      stop("Error reading file: ", e$message)
    }
  )

  # Convert to numeric matrix
  mat <- as.matrix(df)
  if (!is.numeric(mat)) {
    stop("Data could not be coerced to a numeric matrix. Check for non-numeric values inside your matrix.")
  }

  return(mat)
}
