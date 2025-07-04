#' Normalize Gene Expression Data
#'
#' Normalize a numeric gene expression matrix using one of the following methods:
#' log2 transformation, z-score standardization or quantile normalization. The
#' function returns a normalized numeric matrix.
#'
#' @param expr_data A numeric matrix of raw counts (genes x samples).
#' @param method Normalization method: "log2", "zscore" or "quantile".
#'
#' @return A normalized numeric matrix.
#'
#' @importFrom preprocessCore normalize.quantiles
#' @export
#'
#' @examples
#' normalized <- normalization(mat, method = "log2")

normalize_data <- function(expr_data, method = c("log2", "zscore", "quantile")) {
  method <- match.arg(method)
  rn <- rownames(expr_data)
  cn <- colnames(expr_data)

  norm_mat <- switch(
    method,
    # +1 to avoid log(0)
    log2 = log2(expr_data + 1),

    zscore = t(scale(t(expr_data))),

    quantile = normalize.quantiles(expr_data),
  )

  # Restore row and column names
  rownames(norm_mat) <- rn
  colnames(norm_mat) <- cn

  return(norm_mat)
}
