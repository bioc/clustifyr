context("compare_list")
# use capture.output to quiet progress bar from fgsea
shush <- function(...) {
    capture.output(..., file = nullfile())
}

test_that("warning if matrix is not binarized", {
    pbmc_mm <- matrixize_markers(pbmc_markers)
    pbmc_avg <- average_clusters(pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified"
    )
    pbmc_avgb <- binarize_expr(pbmc_avg)
    gene_list_methods <- c("hyper")
    expect_warning(shush(results <- lapply(
        gene_list_methods,
        function(x) {
            compare_lists(pbmc_avg,
                pbmc_mm,
                metric = x
            )
        }
    )))
})


test_that("run all gene list functions", {
    pbmc_mm <- matrixize_markers(pbmc_markers)
    pbmc_avg <- average_clusters(pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified"
    )
    pbmc_avgb <- binarize_expr(pbmc_avg)
    gene_list_methods <- c("spearman", "hyper", "jaccard", "gsea")
    shush(results <- lapply(
        gene_list_methods,
        function(x) {
            compare_lists(pbmc_avgb,
                pbmc_mm,
                metric = x
            )
        }
    ))

    expect_equal(4, length(results))
})

test_that("output intersected genes with details_out option with hyper/jaccard", {
    pbmc_mm <- matrixize_markers(pbmc_markers)
    pbmc_avg <- average_clusters(pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified"
    )
    pbmc_avgb <- binarize_expr(pbmc_avg)
    gene_list_methods <- c("hyper", "jaccard")
    shush(results <- lapply(
        gene_list_methods,
        function(x) {
            compare_lists(pbmc_avgb,
                pbmc_mm,
                metric = x,
                details_out = TRUE
            )
        }
    ))

    expect_equal(2, length(results))
})

test_that("gene list function options", {
    pbmc_mm <- matrixize_markers(pbmc_markers)
    pbmc_avg <- average_clusters(pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified"
    )
    pbmc_avgb <- binarize_expr(pbmc_avg)
    expect_error(suppressWarnings(
        res <- compare_lists(
            pbmc_avgb,
            pbmc_mm,
            metric = "hyper",
            output_high = FALSE,
            n = 5
        )
    ))
})

test_that("run all gene list functions in clustify_lists", {
    gene_list_methods <- c("spearman", "hyper", "jaccard", "gsea")
    shush(results <- lapply(
        gene_list_methods,
        function(x) {
            clustify_lists(
                pbmc_matrix_small,
                per_cell = FALSE,
                metadata = pbmc_meta,
                cluster_col = "classified",
                marker = pbmc_markers,
                marker_inmatrix = FALSE,
                metric = x
            )
        }
    ))

    expect_equal(4, length(results))
})

test_that("gsea outputs in cor matrix format", {
    shush(res <- clustify_lists(
        pbmc_matrix_small,
        per_cell = FALSE,
        metadata = pbmc_meta,
        cluster_col = "classified",
        marker = pbmc_markers,
        marker_inmatrix = FALSE,
        metric = "gsea"
    ))

    res2 <- cor_to_call(res)

    expect_equal(9, nrow(res2))
})

so <- so_pbmc()
test_that("seurat3 object clustify_lists-ing", {
    res <- clustify_lists(
        so,
        per_cell = FALSE,
        marker = pbmc_markers,
        marker_inmatrix = FALSE,
        metric = "jaccard",
        cluster_col = "seurat_clusters",
        seurat_out = TRUE,
        dr = "tsne"
    )
    res <- clustify_lists(
        so,
        per_cell = FALSE,
        marker = pbmc_markers,
        marker_inmatrix = FALSE,
        metric = "jaccard",
        cluster_col = "seurat_clusters",
        seurat_out = FALSE,
        dr = "tsne"
    )
    g <- plot_best_call(
        res,
        seurat_meta(so,
            dr = "tsne"
        ),
        cluster_col = "seurat_clusters",
        plot_r = TRUE,
        x = "tSNE_1",
        y = "tSNE_2"
    )
    expect_true(ggplot2::is.ggplot(g[[1]]))
})

test_that("clustify_lists inserts seurat3 metadata correctly", {
    res <- clustify_lists(
        so,
        per_cell = FALSE,
        marker = pbmc_markers,
        marker_inmatrix = FALSE,
        metric = "jaccard",
        cluster_col = "seurat_clusters",
        seurat_out = TRUE,
        dr = "tsne"
    )
    res2 <- clustify_lists(
        so,
        per_cell = TRUE,
        marker = pbmc_markers,
        marker_inmatrix = FALSE,
        metric = "jaccard",
        cluster_col = "seurat_clusters",
        seurat_out = TRUE,
        dr = "tsne"
    )
    if ("SeuratObject" %in% loadedNamespaces()) {
        expect_true(class(res) %in% c("Seurat"))
    } else {
        expect_true(is.matrix(res))
    }
})

test_that("run all gene list functions and then use consensus_call", {
    pbmc_mm <- matrixize_markers(pbmc_markers)
    pbmc_avg <- average_clusters(pbmc_matrix_small,
        pbmc_meta,
        cluster_col = "classified"
    )
    pbmc_avgb <- binarize_expr(pbmc_avg)
    gene_list_methods <- c("spearman", "hyper", "jaccard", "gsea")
    shush(results <- lapply(
        gene_list_methods,
        function(x) {
            compare_lists(pbmc_avgb,
                pbmc_mm,
                metric = x
            )
        }
    ))
    call_list <- lapply(
        results,
        cor_to_call_rank
    )
    calls <- call_consensus(call_list)
    expect_equal(4, length(results))
})

test_that("run all gene list functions in clustify_lists", {
    res <- clustify_lists(
        pbmc_matrix_small,
        cbmc_m,
        metadata = pbmc_meta,
        cluster_col = "classified",
        metric = "consensus"
    )
    expect_equal(9, nrow(res))
})

test_that("run all gene list functions in clustify_lists and seurat object", {
    res <- clustify_lists(
        so,
        marker = cbmc_m,
        dr = "tsne",
        cluster_col = "seurat_clusters",
        metric = "consensus",
        seurat_out = TRUE
    )
    expect_true(is.data.frame(res) | "Seurat" %in% class(res))
})

test_that("lists of genes will work with posneg", {
    lst_of_markers <- split(pbmc_markers$gene, pbmc_markers$cluster)
    res <- clustify_lists(
        input = pbmc_matrix_small,
        per_cell = FALSE,
        cluster_col = "classified",
        metadata = pbmc_meta,
        marker = lst_of_markers,
        marker_inmatrix = TRUE,
        metric = "posneg",
        seurat_out = FALSE
    )
    expect_true(ncol(res) == length(lst_of_markers))
})

test_that("clustify_lists input_markers mode", {
    pbmc_mm <- matrixize_markers(pbmc_markers)
    pbmc_input_mm <- pos_neg_marker(pbmc_mm[1:3, ])
    results <- lapply(
        c("hyper", "spearman"),
        function(x) {
            clustify_lists(pbmc_input_mm,
                pbmc_mm,
                metric = x,
                input_markers = TRUE
            )
        }
    )
    expect_equal(2, length(results))
})

test_that("clustify_lists input_markers mode with uneven number of marker per cluster", {
    pbmc_mm <- matrixize_markers(pbmc_markers)
    pbmc_input_mm <- pos_neg_marker(pbmc_mm[1:3, ])
    results <- lapply(
        c("jaccard"),
        function(x) {
            clustify_lists(pbmc_input_mm,
                split(pbmc_markers$gene, pbmc_markers$cluster),
                metric = x,
                input_markers = TRUE
            )
        }
    )
    expect_equal(1, length(results))
})
