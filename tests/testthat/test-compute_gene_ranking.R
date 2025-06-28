# Correct output in tems of format and sorting
# Gene_cluster is not named
# Gene_clusters names not in expr_mat rownames -> test with gene X
# No genes in the cluster
# Maybe check rownames, TBD


test_that("compute_gene_ranking returns correct output format and sorting", {
  expr_mat <- matrix(rnorm(30), nrow = 6, ncol = 5)
  rownames(expr_mat) <- paste0("Gene", 1:6)

  gene_clusters <- c(1,1,1,2,2,2)
  names(gene_clusters) <- rownames(expr_mat)

  ranking <- compute_gene_ranking(expr_mat, gene_clusters, 1, 2)

  expect_true(is.numeric(ranking))
  expect_true(!is.null(names(ranking)))
  expect_true(all(names(ranking) %in% rownames(expr_mat)))
  expect_equal(sort(ranking, decreasing = TRUE), ranking)
})

test_that("error if gene_clusters is not named", {
  expr_mat <- matrix(rnorm(10), nrow = 2)
  rownames(expr_mat) <- c("Gene1", "Gene2")

  gene_clusters <- c(1, 2)

  expect_error(compute_gene_ranking(expr_mat, gene_clusters, 1, 2),
               "gene_clusters must be a named vector")
})

test_that("error if gene_clusters names not in expr_mat rownames", {
  expr_mat <- matrix(rnorm(10), nrow = 2)
  rownames(expr_mat) <- c("Gene1", "Gene2")

  gene_clusters <- c(1, 2)
  names(gene_clusters) <- c("Gene1", "GeneX")  # "GeneX" does not exist

  expect_error(compute_gene_ranking(expr_mat, gene_clusters, 1, 2),
               "All gene_clusters names must be in expr_mat rownames")
})

test_that("works when one cluster has no genes (returns numeric vector)", {
  expr_mat <- matrix(rnorm(15), nrow = 5)
  rownames(expr_mat) <- paste0("Gene", 1:5)

  gene_clusters <- c(1, 1, 1, 1, 1)
  names(gene_clusters) <- rownames(expr_mat)

  ranking <- compute_gene_ranking(expr_mat, gene_clusters, 1, 2)

  expect_true(is.numeric(ranking))
})

