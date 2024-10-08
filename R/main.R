#' Compare scRNA-seq data to reference data.
#'
#' @export
clustify <- function(input, ...) {
    UseMethod("clustify", input)
}

#' @rdname clustify
#' @param input single-cell expression matrix or Seurat object
#' @param metadata cell cluster assignments,
#'   supplied as a vector or data.frame.
#'   If data.frame is supplied then `cluster_col` needs to be set.
#'   Not required if running correlation per cell.
#' @param ref_mat reference expression matrix
#' @param cluster_col column in metadata that contains cluster ids per cell.
#'   Will default to first column of metadata if not supplied.
#'   Not required if running correlation per cell.
#' @param query_genes A vector of genes of interest to compare. If NULL, then
#'   common genes between the expr_mat and ref_mat
#'   will be used for comparision.
#' @param n_genes number of genes limit for Seurat variable genes, by default 1000,
#'   set to 0 to use all variable genes (generally not recommended)
#' @param per_cell if true run per cell, otherwise per cluster.
#' @param n_perm number of permutations, set to 0 by default
#' @param compute_method method(s) for computing similarity scores
#' @param pseudobulk_method method used for summarizing clusters, options are mean (default), median, truncate (10% truncated mean), or trimean, max, min
#' @param use_var_genes if providing a seurat object, use the variable genes
#'   (stored in seurat_object@var.genes) as the query_genes.
#' @param dr stored dimension reduction
#' @param seurat_out output cor matrix or called seurat object
#'  (deprecated, use obj_out instead)
#' @param verbose whether to report certain variables chosen and steps
#' @param lookuptable if not supplied, will look in built-in table
#'  for object parsing
#' @param rm0 consider 0 as missing data, recommended for per_cell
#' @param obj_out whether to output object instead of cor matrix
#' @param vec_out only output a result vector in the same order as metadata
#' @param rename_prefix prefix to add to type and r column names
#' @param threshold identity calling minimum correlation score threshold,
#'  only used when obj_out = TRUE
#' @param low_threshold_cell option to remove clusters with too few cells
#' @param exclude_genes a vector of gene names to throw out of query
#' @param if_log input data is natural log,
#' averaging will be done on unlogged data
#' @param organism for GO term analysis, organism name: human - 'hsapiens', mouse - 'mmusculus'
#' @param plot_name name for saved pdf, if NULL then no file is written (default)
#' @param rds_name name for saved rds of rank_diff, if NULL then no file is written (default)
#' @param expand_unassigned test all ref clusters for unassigned results
#' @param ... additional arguments to pass to compute_method function
#'
#' @return single cell object with identity assigned in metadata,
#'   or matrix of correlation values, clusters from input as row names, cell
#'   types from ref_mat as column names
#'
#' @examples
#' # Annotate a matrix and metadata
#' clustify(
#'     input = pbmc_matrix_small,
#'     metadata = pbmc_meta,
#'     ref_mat = cbmc_ref,
#'     query_genes = pbmc_vargenes,
#'     cluster_col = "RNA_snn_res.0.5",
#'     verbose = TRUE
#' )
#'
#' # Annotate using a different method
#' clustify(
#'     input = pbmc_matrix_small,
#'     metadata = pbmc_meta,
#'     ref_mat = cbmc_ref,
#'     query_genes = pbmc_vargenes,
#'     cluster_col = "RNA_snn_res.0.5",
#'     compute_method = "cosine"
#' )
#'
#' # Annotate a SingleCellExperiment object
#' sce <- sce_pbmc()
#' clustify(
#'     sce,
#'     cbmc_ref,
#'     cluster_col = "clusters",
#'     obj_out = TRUE,
#'     per_cell = FALSE,
#'     dr = "umap"
#' )
#'
#' # Annotate a Seurat object
#' so <- so_pbmc()
#' clustify(
#'     so,
#'     cbmc_ref,
#'     cluster_col = "seurat_clusters",
#'     obj_out = TRUE,
#'     per_cell = FALSE,
#'     dr = "umap"
#' )
#'
#' # Annotate (and return) a Seurat object per-cell
#' clustify(
#'     input = so,
#'     ref_mat = cbmc_ref,
#'     cluster_col = "seurat_clusters",
#'     obj_out = TRUE,
#'     per_cell = TRUE,
#'     dr = "umap"
#' )
#' @export
clustify.default <- function(
    input,
    ref_mat,
    metadata = NULL,
    cluster_col = NULL,
    query_genes = NULL,
    n_genes = 1000,
    per_cell = FALSE,
    n_perm = 0,
    compute_method = "spearman",
    pseudobulk_method = "mean",
    verbose = TRUE,
    lookuptable = NULL,
    rm0 = FALSE,
    obj_out = TRUE,
    seurat_out = obj_out,
    vec_out = FALSE,
    rename_prefix = NULL,
    threshold = "auto",
    low_threshold_cell = 0,
    exclude_genes = c(),
    if_log = TRUE,
    organism = "hsapiens",
    plot_name = NULL,
    rds_name = NULL,
    expand_unassigned = FALSE,
    ...) {
    if (!compute_method %in% clustifyr_methods) {
        stop(compute_method, " correlation method not implemented",
            call. = FALSE
        )
    }

    input_original <- input
    if (!inherits(input_original, c("matrix", "Matrix", "data.frame"))) {
        temp <- parse_loc_object(
            input,
            type = class(input),
            expr_loc = NULL,
            meta_loc = NULL,
            var_loc = NULL,
            cluster_col = cluster_col,
            lookuptable = lookuptable
        )

        if (!(is.null(temp[["expr"]]))) {
            message("recognized object type - ", class(input))
        }

        input <- temp[["expr"]]
        metadata <- temp[["meta"]]
        if (is.null(query_genes)) {
            query_genes <- temp[["var"]]
        }

        if (is.null(cluster_col)) {
            cluster_col <- temp[["col"]]
        }
    }

    if (is.null(metadata) && !per_cell) {
        stop("`metadata` needed for per cluster analysis", call. = FALSE)
    }

    if (!is.null(cluster_col) &&
        !cluster_col %in% colnames(metadata)) {
        stop("given `cluster_col` is not a column in `metadata`", call. = FALSE)
    }

    if (is.null(query_genes) || length(query_genes) == 0) {
        message(
            "Variable features not available, using all genes instead\n",
            "consider supplying variable features to `query_genes` argument."
        )
        query_genes <- NULL
    }

    expr_mat <- input

    # select gene subsets
    gene_constraints <- get_common_elements(
        rownames(expr_mat),
        rownames(ref_mat),
        query_genes
    )

    if (length(exclude_genes) > 0) {
        gene_constraints <- setdiff(gene_constraints, exclude_genes)
    }

    if (verbose) {
        message("using # of genes: ", length(gene_constraints))
        if (length(gene_constraints) >= 2000) {
            message(
                "using a high number (>=2000) genes to calculate correlation ",
                "please consider feature selection to improve performance"
            )
        }
    }

    expr_mat <- expr_mat[gene_constraints, , drop = FALSE]
    ref_mat <- ref_mat[gene_constraints, , drop = FALSE]

    if (!per_cell) {
        if (is.vector(metadata)) {
            cluster_ids <- metadata
        } else if (is.factor(metadata)) {
            cluster_ids <- as.character(metadata)
        } else if (is.data.frame(metadata) & !is.null(cluster_col)) {
            cluster_ids <- metadata[[cluster_col]]
        } else {
            stop(
                "metadata not formatted correctly,
           supply either a character vector or a dataframe",
                call. = FALSE
            )
        }
        if (is.factor(cluster_ids)) {
            cluster_ids <- as.character(cluster_ids)
        }
        cluster_ids[is.na(cluster_ids)] <- "orig.NA"
    }

    if (per_cell) {
        cluster_ids <- colnames(expr_mat)
    }

    if (n_perm == 0) {
        res <- get_similarity(
            expr_mat,
            ref_mat,
            cluster_ids = cluster_ids,
            per_cell = per_cell,
            compute_method = compute_method,
            pseudobulk_method = pseudobulk_method,
            rm0 = rm0,
            if_log = if_log,
            low_threshold = low_threshold_cell,
            ...
        )
    } else {
        # run permutation
        res <- permute_similarity(
            expr_mat,
            ref_mat,
            cluster_ids = cluster_ids,
            n_perm = n_perm,
            per_cell = per_cell,
            compute_method = compute_method,
            pseudobulk_method = pseudobulk_method,
            rm0 = rm0,
            ...
        )
    }

    if (verbose) {
        message("similarity computation completed, matrix of ", dim(res)[1], " x ", dim(res)[2], ", preparing output")
    }

    obj_out <- seurat_out
    if (obj_out &&
        !inherits(input_original, c(
            "matrix",
            "Matrix",
            "data.frame"
        )) || (vec_out &&
        inherits(input_original, c(
            "matrix",
            "Matrix",
            "data.frame"
        )))) {
        df_temp <- cor_to_call(
            res,
            metadata = metadata,
            cluster_col = cluster_col,
            threshold = threshold
        )

        df_temp_full <- call_to_metadata(
            df_temp,
            metadata = metadata,
            cluster_col = cluster_col,
            per_cell = per_cell,
            rename_prefix = rename_prefix
        )

        if (vec_out) {
            if (is.null(rename_prefix)) {
                return(df_temp_full[["type"]])
            } else {
                return(df_temp_full[[paste0(rename_prefix, "_type")]])
            }
        }

        out <- insert_meta_object(input_original,
            df_temp_full,
            lookuptable = lookuptable
        )

        if (!is.null(plot_name)) {
            message("saving rank diff plot")
            avg_mat <- average_clusters(
                expr_mat,
                cluster_ids
            )
            assess_rank_bias(
                avg_mat = avg_mat,
                ref_mat = ref_mat,
                query_genes = query_genes,
                res = df_temp,
                organism = organism,
                plot_name = plot_name,
                rds_name = rds_name,
                expand_unassigned = expand_unassigned
            )
        }
        return(out)
    } else {
        return(res)
    }
}

#' @rdname clustify
#' @export
clustify.Seurat <- function(
    input,
    ref_mat,
    cluster_col = NULL,
    query_genes = NULL,
    n_genes = 1000,
    per_cell = FALSE,
    n_perm = 0,
    compute_method = "spearman",
    pseudobulk_method = "mean",
    use_var_genes = TRUE,
    dr = "umap",
    obj_out = TRUE,
    seurat_out = obj_out,
    vec_out = FALSE,
    threshold = "auto",
    verbose = TRUE,
    rm0 = FALSE,
    rename_prefix = NULL,
    exclude_genes = c(),
    metadata = NULL,
    organism = "hsapiens",
    plot_name = NULL,
    rds_name = NULL,
    expand_unassigned = FALSE,
    ...) {
    s_object <- input
    # for seurat 3.0 +
    expr_mat <- object_data(s_object, "data")
    vec <- FALSE
    if (!is.null(metadata)) {
        if (is.vector(metadata)) {
            vec <- TRUE
        } else if (is.factor(metadata)) {
            vec <- TRUE
            metadata <- as.character(metadata)
        }
    } else {
        metadata <- seurat_meta(s_object, dr = dr)
    }

    if (use_var_genes && is.null(query_genes)) {
        query_genes <- object_data(s_object, "var.genes", n_genes)
    }

    if (verbose) {
        message("object data retrieval complete, moving to similarity computation")
    }

    res <- clustify(
        expr_mat,
        ref_mat,
        metadata,
        query_genes,
        per_cell = per_cell,
        n_perm = n_perm,
        cluster_col = cluster_col,
        compute_method = compute_method,
        pseudobulk_method = pseudobulk_method,
        verbose = verbose,
        rm0 = rm0,
        exclude_genes = exclude_genes,
        ...
    )

    if (n_perm != 0) {
        res <- -log(res$p_val + .01, 10)
    }
    obj_out <- seurat_out
    if (!obj_out && !vec_out || vec) {
        res
    } else {
        df_temp <- cor_to_call(
            res,
            metadata = metadata,
            cluster_col = cluster_col,
            threshold = threshold
        )

        df_temp_full <- call_to_metadata(
            df_temp,
            metadata = metadata,
            cluster_col = cluster_col,
            per_cell = per_cell,
            rename_prefix = rename_prefix
        )

        if (!is.null(plot_name)) {
            message("saving rank diff plot")
            avg_mat <- average_clusters(
                expr_mat,
                metadata,
                cluster_col
            )
            assess_rank_bias(
                avg_mat = avg_mat,
                ref_mat = ref_mat,
                query_genes = query_genes,
                res = df_temp,
                organism = organism,
                plot_name = plot_name,
                rds_name = rds_name,
                expand_unassigned = expand_unassigned
            )
        }

        if (vec_out) {
            if (is.null(rename_prefix)) {
                return(df_temp_full[["type"]])
            } else {
                return(df_temp_full[[paste0(rename_prefix, "_type")]])
            }
        }

        if ("SeuratObject" %in% loadedNamespaces()) {
            s_object <- write_meta(s_object, df_temp_full)
            return(s_object)
        } else {
            message("seurat not loaded, returning cor_mat instead")
            return(res)
        }
        s_object
    }
}

#' @rdname clustify
#' @export
clustify.SingleCellExperiment <- function(
    input,
    ref_mat,
    cluster_col = NULL,
    query_genes = NULL,
    per_cell = FALSE,
    n_perm = 0,
    compute_method = "spearman",
    pseudobulk_method = "mean",
    use_var_genes = TRUE,
    dr = "umap",
    obj_out = TRUE,
    seurat_out = obj_out,
    vec_out = FALSE,
    threshold = "auto",
    verbose = TRUE,
    rm0 = FALSE,
    rename_prefix = NULL,
    exclude_genes = c(),
    metadata = NULL,
    organism = "hsapiens",
    plot_name = NULL,
    rds_name = NULL,
    expand_unassigned = FALSE,
    ...) {
    s_object <- input
    expr_mat <- object_data(s_object, "data")
    vec <- FALSE
    if (!is.null(metadata)) {
        if (is.vector(metadata)) {
            vec <- TRUE
        } else if (is.factor(metadata)) {
            vec <- TRUE
            metadata <- as.character(metadata)
        }
    } else {
        metadata <- object_data(s_object, "meta.data")
    }


    if (verbose) {
        message("object data retrieval complete, moving to similarity computation")
    }

    res <- clustify(
        expr_mat,
        ref_mat,
        metadata,
        query_genes,
        per_cell = per_cell,
        n_perm = n_perm,
        cluster_col = cluster_col,
        compute_method = compute_method,
        pseudobulk_method = pseudobulk_method,
        verbose = verbose,
        rm0 = rm0,
        exclude_genes = exclude_genes,
        ...
    )

    if (n_perm != 0) {
        res <- -log(res$p_val + .01, 10)
    }
    obj_out <- seurat_out
    if (!obj_out && !vec_out) {
        res
    } else {
        df_temp <- cor_to_call(
            res,
            metadata = metadata,
            cluster_col = cluster_col,
            threshold = threshold
        )

        df_temp_full <- call_to_metadata(
            df_temp,
            metadata = metadata,
            cluster_col = cluster_col,
            per_cell = per_cell,
            rename_prefix = rename_prefix
        )

        if (!is.null(plot_name)) {
            message("saving rank diff plot")
            avg_mat <- average_clusters(
                expr_mat,
                metadata,
                cluster_col
            )
            assess_rank_bias(
                avg_mat = avg_mat,
                ref_mat = ref_mat,
                query_genes = query_genes,
                res = df_temp,
                organism = organism,
                plot_name = plot_name,
                rds_name = rds_name,
                expand_unassigned = expand_unassigned
            )
        }

        if (vec_out) {
            if (is.null(rename_prefix)) {
                return(df_temp_full[["type"]])
            } else {
                return(df_temp_full[[paste0(rename_prefix, "_type")]])
            }
        }

        if ("SingleCellExperiment" %in% loadedNamespaces()) {
            if (!(is.null(rename_prefix))) {
                col_type <- stringr::str_c(rename_prefix, "_type")
                col_r <- stringr::str_c(rename_prefix, "_r")
            } else {
                col_type <- "type"
                col_r <- "r"
            }
            colDatatemp <- metadata
            colDatatemp[[col_type]] <- df_temp_full[[col_type]]
            colDatatemp[[col_r]] <- df_temp_full[[col_r]]
            s_object <- write_meta(s_object, colDatatemp)
            return(s_object)
        } else {
            message("SingleCellExperiment not loaded, returning cor_mat instead")
            return(res)
        }
        s_object
    }
}

#' Correlation functions available in clustifyr
#' @examples
#' clustifyr_methods
#' @export
clustifyr_methods <- c(
    "pearson",
    "spearman",
    "cosine",
    "kl_divergence",
    "kendall"
)

#' Main function to compare scRNA-seq data to gene lists.
#'
#' @export
clustify_lists <- function(input, ...) {
    UseMethod("clustify_lists", input)
}

#' @rdname clustify_lists
#' @param input single-cell expression matrix, Seurat object, or SingleCellExperiment
#' @param marker matrix or dataframe of candidate genes for each cluster
#' @param marker_inmatrix whether markers genes are already in preprocessed
#'   matrix form
#' @param metadata cell cluster assignments,
#'   supplied as a vector or data.frame.
#'   If data.frame is supplied then `cluster_col` needs to be set.
#'   Not required if running correlation per cell.
#' @param cluster_col column in metadata with cluster number
#' @param if_log input data is natural log, averaging will be done on
#'   unlogged data
#' @param per_cell compare per cell or per cluster
#' @param topn number of top expressing genes to keep from input matrix
#' @param cut expression cut off from input matrix
#' @param genome_n number of genes in the genome
#' @param metric adjusted p-value for hypergeometric test, or jaccard index
#' @param output_high if true (by default to fit with rest of package),
#' -log10 transform p-value
#' @param lookuptable if not supplied, will look in built-in table
#'  for object parsing
#' @param obj_out whether to output object instead of cor matrix
#' @param vec_out only output a result vector in the same order as metadata
#' @param rename_prefix prefix to add to type and r column names
#' @param threshold identity calling minimum correlation score threshold,
#' only used when obj_out = T
#' @param low_threshold_cell option to remove clusters with too few cells
#' @param dr stored dimension reduction
#' @param seurat_out output cor matrix or called seurat object
#'   (deprecated, use obj_out instead)
#' @param verbose whether to report certain variables chosen and steps
#' @param input_markers whether input is marker data.frame of 0 and 1s (output of pos_neg_marker), and uses alternate enrichment mode
#' @param details_out whether to also output shared gene list from jaccard
#' @param ... passed to matrixize_markers
#' @examples
#' # Annotate a matrix and metadata
# clustify_lists(
#     input = pbmc_matrix_small,
#     marker = cbmc_m,
#     metadata = pbmc_meta,
#     cluster_col = "classified",
#     verbose = TRUE
# )
#'
#' # Annotate using a different method
#' clustify_lists(
#'     input = pbmc_matrix_small,
#'     marker = cbmc_m,
#'     metadata = pbmc_meta,
#'     cluster_col = "classified",
#'     verbose = TRUE,
#'     metric = "jaccard"
#' )
#' @return matrix of numeric values, clusters from input as row names,
#' cell types from marker_mat as column names

#' @export
clustify_lists.default <- function(
    input,
    marker,
    marker_inmatrix = TRUE,
    metadata = NULL,
    cluster_col = NULL,
    if_log = TRUE,
    per_cell = FALSE,
    topn = 800,
    cut = 0,
    genome_n = 30000,
    metric = "hyper",
    output_high = TRUE,
    lookuptable = NULL,
    obj_out = TRUE,
    seurat_out = obj_out,
    vec_out = FALSE,
    rename_prefix = NULL,
    threshold = 0,
    low_threshold_cell = 0,
    verbose = TRUE,
    input_markers = FALSE,
    details_out = FALSE,
    ...) {
    input_original <- input
    if (!inherits(input, c("matrix", "Matrix", "data.frame"))) {
        temp <- parse_loc_object(
            input,
            type = class(input),
            expr_loc = NULL,
            meta_loc = NULL,
            var_loc = NULL,
            cluster_col = cluster_col,
            lookuptable = lookuptable
        )
        input <- temp[["expr"]]
        metadata <- temp[["meta"]]
        cluster_info <- metadata
        if (is.null(cluster_col)) {
            cluster_col <- temp[["col"]]
        }
    } else {
        cluster_info <- metadata
    }

    if (metric %in% c("posneg", "pct")) {
        per_cell <- TRUE
    }
    if (input_markers) {
        per_cell <- TRUE
    }
    if (!(per_cell)) {
        input <- average_clusters(input,
            cluster_info,
            if_log = if_log,
            cluster_col = cluster_col,
            low_threshold = low_threshold_cell
        )
    }

    if (!input_markers) {
        bin_input <- binarize_expr(input, n = topn, cut = cut)
    } else {
        bin_input <- as.matrix(input_original)
    }

    if (marker_inmatrix != TRUE & metric != "posneg") {
        marker <- matrixize_markers(
            marker,
            ...
        )
        if (verbose) {
            message("number of total markers: ", nrow(marker))
        }
    }

    if (metric == "consensus") {
        results <- lapply(
            c("hyper", "jaccard", "pct", "posneg"),
            function(x) {
                clustify_lists(
                    input_original,
                    marker,
                    metadata = cluster_info,
                    cluster_col = cluster_col,
                    metric = x
                )
            }
        )
        call_list <- lapply(
            results,
            cor_to_call_rank
        )
        res <- call_consensus(call_list)
    } else if (metric == "pct") {
        res <- gene_pct_markerm(input,
            marker,
            cluster_info,
            cluster_col = cluster_col
        )
    } else if (metric == "gsea") {
        res <- compare_lists(
            input,
            marker_mat = marker,
            n = genome_n,
            metric = "gsea",
            output_high = output_high
        )
    } else if (metric != "posneg") {
        res <- compare_lists(
            bin_input,
            marker_mat = marker,
            n = genome_n,
            metric = metric,
            output_high = output_high,
            details_out = details_out
        )
    } else {
        if (is.data.frame(marker)) {
            marker <- as.matrix(marker)
        }
        if (!is.numeric(marker)) {
            marker <- pos_neg_marker(marker)
        }
        res <- pos_neg_select(input,
            marker,
            cluster_info,
            cluster_col = cluster_col,
            cutoff_score = NULL
        )
    }

    if (verbose) {
        message("similarity computation completed, matrix of ", dim(res)[1], " x ", dim(res)[2], ", preparing output")
    }
    obj_out <- seurat_out
    if ((!inherits(input_original, c("matrix", "Matrix", "data.frame")) &&
        obj_out) || (vec_out &&
        inherits(input_original, c(
            "matrix",
            "Matrix",
            "data.frame"
        )))) {
        if (metric != "consensus") {
            df_temp <- cor_to_call(
                res,
                metadata = metadata,
                cluster_col = cluster_col,
                threshold = threshold
            )

            df_temp_full <- call_to_metadata(
                df_temp,
                metadata = metadata,
                cluster_col = cluster_col,
                per_cell = per_cell,
                rename_prefix = rename_prefix
            )
        } else {
            df_temp_full <- res
        }

        if (vec_out) {
            if (is.null(rename_prefix)) {
                return(df_temp_full[["type"]])
            } else {
                return(df_temp_full[[paste0(rename_prefix, "_type")]])
            }
        }

        out <- insert_meta_object(input_original,
            df_temp_full,
            lookuptable = lookuptable
        )

        return(out)
    } else {
        return(res)
    }
}

#' @rdname clustify_lists
#' @export
clustify_lists.Seurat <- function(
    input,
    metadata = NULL,
    cluster_col = NULL,
    if_log = TRUE,
    per_cell = FALSE,
    topn = 800,
    cut = 0,
    marker,
    marker_inmatrix = TRUE,
    genome_n = 30000,
    metric = "hyper",
    output_high = TRUE,
    dr = "umap",
    obj_out = TRUE,
    seurat_out = obj_out,
    vec_out = FALSE,
    threshold = 0,
    rename_prefix = NULL,
    verbose = TRUE,
    details_out = FALSE,
    ...) {
    s_object <- input
    # for seurat 3.0 +
    input <- object_data(s_object, "data")
    vec <- FALSE
    if (!is.null(metadata)) {
        if (is.vector(metadata)) {
            vec <- TRUE
        } else if (is.factor(metadata)) {
            vec <- TRUE
            metadata <- as.character(metadata)
        }
    } else {
        metadata <- object_data(s_object, "meta.data")
    }
    cluster_info <- metadata

    if (verbose) {
        message("object data retrieval complete, moving to similarity computation")
    }

    res <- clustify_lists(
        input,
        per_cell = per_cell,
        metadata = cluster_info,
        if_log = if_log,
        cluster_col = cluster_col,
        topn = topn,
        cut = cut,
        marker,
        marker_inmatrix = marker_inmatrix,
        genome_n = genome_n,
        metric = metric,
        output_high = output_high,
        verbose = verbose,
        details_out = details_out,
        ...
    )
    obj_out <- seurat_out
    if (!obj_out && !vec_out || vec) {
        res
    } else {
        if (metric != "consensus") {
            df_temp <- cor_to_call(
                res,
                metadata = metadata,
                cluster_col = cluster_col,
                threshold = threshold
            )
        } else {
            df_temp <- res
            colnames(df_temp)[1] <- cluster_col
        }

        df_temp_full <- call_to_metadata(
            df_temp,
            metadata = metadata,
            cluster_col = cluster_col,
            per_cell = per_cell,
            rename_prefix = rename_prefix
        )

        if (vec_out) {
            if (is.null(rename_prefix)) {
                return(df_temp_full[["type"]])
            } else {
                return(df_temp_full[[paste0(rename_prefix, "_type")]])
            }
        }

        if ("SeuratObject" %in% loadedNamespaces()) {
            s_object <- write_meta(s_object, df_temp_full)
            return(s_object)
        } else {
            message("seurat not loaded, returning cor_mat instead")
            return(res)
        }
        s_object
    }
}

#' @rdname clustify_lists
#' @export
clustify_lists.SingleCellExperiment <- function(
    input,
    metadata = NULL,
    cluster_col = NULL,
    if_log = TRUE,
    per_cell = FALSE,
    topn = 800,
    cut = 0,
    marker,
    marker_inmatrix = TRUE,
    genome_n = 30000,
    metric = "hyper",
    output_high = TRUE,
    dr = "umap",
    obj_out = TRUE,
    seurat_out = obj_out,
    vec_out = FALSE,
    threshold = 0,
    rename_prefix = NULL,
    verbose = TRUE,
    details_out = FALSE,
    ...) {
    s_object <- input
    expr_mat <- object_data(s_object, "data")
    vec <- FALSE
    if (!is.null(metadata)) {
        if (is.vector(metadata)) {
            vec <- TRUE
        } else if (is.factor(metadata)) {
            vec <- TRUE
            metadata <- as.character(metadata)
        }
    } else {
        metadata <- object_data(s_object, "meta.data")
    }

    if (verbose) {
        message("object data retrieval complete, moving to similarity computation")
    }

    res <- clustify_lists(
        expr_mat,
        per_cell = per_cell,
        metadata = metadata,
        if_log = if_log,
        cluster_col = cluster_col,
        topn = topn,
        cut = cut,
        marker,
        marker_inmatrix = marker_inmatrix,
        genome_n = genome_n,
        metric = metric,
        output_high = output_high,
        details_out = details_out,
        ...
    )

    if (!obj_out && !vec_out || vec) {
        res
    } else {
        df_temp <- cor_to_call(
            res,
            metadata = metadata,
            cluster_col = cluster_col,
            threshold = threshold
        )

        df_temp_full <- call_to_metadata(
            df_temp,
            metadata = metadata,
            cluster_col = cluster_col,
            per_cell = per_cell,
            rename_prefix = rename_prefix
        )

        if (vec_out) {
            if (is.null(rename_prefix)) {
                return(df_temp_full[["type"]])
            } else {
                return(df_temp_full[[paste0(rename_prefix, "_type")]])
            }
        }

        if ("SingleCellExperiment" %in% loadedNamespaces()) {
            if (!(is.null(rename_prefix))) {
                col_type <- stringr::str_c(rename_prefix, "_type")
                col_r <- stringr::str_c(rename_prefix, "_r")
            } else {
                col_type <- "type"
                col_r <- "r"
            }
            colDatatemp <- metadata
            colDatatemp[[col_type]] <- df_temp_full[[col_type]]
            colDatatemp[[col_r]] <- df_temp_full[[col_r]]
            s_object <- write_meta(s_object, colDatatemp)
            return(s_object)
        } else {
            message("SingleCellExperiment not loaded, returning cor_mat instead")
            return(res)
        }
        s_object
    }
}
