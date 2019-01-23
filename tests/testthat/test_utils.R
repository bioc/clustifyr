context("clustify")

test_that("get_vargenes works for both matrix and dataframe form", {
  pbmc4k_mm <- matrixize_markers(pbmc4k_markers)
  var1 <- get_vargenes(pbmc4k_mm)
  var2 <- get_vargenes(pbmc4k_markers)

  expect_equal(var1[1], var2[1])
})

test_that("average_clusters works as intended", {
  pbmc4k_avg2 <- average_clusters(pbmc4k_matrix, pbmc4k_meta)
  expect_equal(nrow(pbmc4k_avg2), nrow(pbmc4k_avg))
})

test_that("average_clusters works when one cluster contains only 1 cell", {
  pbmc4k_meta2 <- pbmc4k_meta
  pbmc4k_meta2[1,"cluster"] <- 15
  pbmc4k_avg2 <- average_clusters(pbmc4k_matrix, pbmc4k_meta2)
  expect_equal(ncol(pbmc4k_avg2), ncol(pbmc4k_avg) + 1)
})

test_that("average_clusters_filter works on strings", {
  avg1 <- average_clusters_filter(pbmc4k_matrix, pbmc4k_meta,
                                  filter_on = "cluster",
                                  filter_method = "==",
                                  filter_value = "1")
  remove_background(pbmc4k_matrix, avg1, 1)
  expect_equal(class(avg1), "numeric")
})

test_that("cor_to_call threshold works as intended", {
  res <- clustify(
    input = pbmc4k_matrix,
    metadata = pbmc4k_meta,
    bulk_mat = pbmc_bulk_matrix,
    query_genes = pbmc4k_vargenes,
    cluster_col = "cluster"
  )
  call1 <- cor_to_call(res,
              metadata = pbmc4k_meta,
              col = "cluster",
              collapse_to_cluster = FALSE,
              threshold = 0.5)

  expect_true("r<0.5, unassigned" %in% call1$type)
})