#' Compute similarity of matrices
#'
#' @param expr_mat single-cell expression matrix
#' @param ref_mat reference expression matrix
#' @param cluster_ids vector of cluster ids for each cell
#' @param compute_method method(s) for computing similarity scores
#' @param pseudobulk_method method used for summarizing clusters, options are mean (default), median, truncate (10% truncated mean), or trimean, max, min
#' @param per_cell run per cell?
#' @param rm0 consider 0 as missing data, recommended for per_cell
#' @param if_log input data is natural log,
#' averaging will be done on unlogged data
#' @param low_threshold option to remove clusters with too few cells
#' @param ... additional parameters not used yet
#' @return matrix of numeric values, clusters from expr_mat as row names,
#'  cell types from ref_mat as column names
get_similarity <- function(
    expr_mat,
    ref_mat,
    cluster_ids,
    compute_method,
    pseudobulk_method = "mean",
    per_cell = FALSE,
    rm0 = FALSE,
    if_log = TRUE,
    low_threshold = 0,
    ...) {
    if (nrow(expr_mat) == 0) {
        stop("after subsetting to shared genes",
            "query expression matrix has 0 rows",
            call. = FALSE
        )
    }

    if (ncol(expr_mat) == 0) {
        stop("query expression matrix has 0 cols",
            call. = FALSE
        )
    }

    if (nrow(ref_mat) == 0) {
        stop("after subsetting to shared genes",
            "reference expression matrix has 0 rows",
            call. = FALSE
        )
    }

    if (ncol(ref_mat) == 0) {
        stop("reference expression matrix has 0 cols",
            call. = FALSE
        )
    }

    ref_clust <- colnames(ref_mat)
    if (ncol(expr_mat) != length(cluster_ids)) {
        stop("number of cells in expression matrix not equal",
            "to metadata/cluster_col",
            call. = FALSE
        )
    }

    if (sum(is.na(cluster_ids)) > 0) {
        message("reassigning NAs to unknown")
        cluster_ids <- factor(cluster_ids)
        cluster_ids <-
            factor(
                cluster_ids,
                levels = c(levels(cluster_ids), NA),
                labels = c(levels(cluster_ids), "unknown"),
                exclude = NULL
            )
        cluster_ids <- as.character(cluster_ids)
    }

    if (!per_cell) {
        sc_clust <- sort(unique(cluster_ids))
        clust_avg <- average_clusters(
            expr_mat,
            cluster_ids,
            if_log = if_log,
            low_threshold = low_threshold,
            method = pseudobulk_method
        )
    } else {
        sc_clust <- cluster_ids
        clust_avg <- expr_mat
    }

    assigned_score <- calc_similarity(
        clust_avg,
        ref_mat,
        compute_method,
        rm0 = rm0,
        ...
    )

    if (low_threshold == 0) {
        rownames(assigned_score) <- sc_clust
        colnames(assigned_score) <- ref_clust
    }

    return(assigned_score)
}

#' Compute a p-value for similarity using permutation
#'
#' @description Permute cluster labels to calculate empirical p-value
#'
#'
#' @param expr_mat single-cell expression matrix
#' @param cluster_ids clustering info of single-cell data assume that
#'  genes have ALREADY BEEN filtered
#' @param ref_mat reference expression matrix
#' @param n_perm number of permutations
#' @param per_cell run per cell?
#' @param compute_method method(s) for computing similarity scores
#' @param pseudobulk_method method used for summarizing clusters, options are mean (default), median, truncate (10% truncated mean), or trimean, max, min
#' @param rm0 consider 0 as missing data, recommended for per_cell
#' @param ... additional parameters
#' @return matrix of numeric values
permute_similarity <- function(
    expr_mat,
    ref_mat,
    cluster_ids,
    n_perm,
    per_cell = FALSE,
    compute_method,
    pseudobulk_method = "mean",
    rm0 = FALSE,
    ...) {
    ref_clust <- colnames(ref_mat)

    if (!per_cell) {
        sc_clust <- sort(unique(cluster_ids))
        clust_avg <- average_clusters(
            expr_mat,
            cluster_ids,
            method = pseudobulk_method
        )
    } else {
        sc_clust <- colnames(expr_mat)
        clust_avg <- expr_mat
    }

    assigned_score <- calc_similarity(clust_avg,
        ref_mat,
        compute_method,
        rm0 = rm0,
        ...
    )

    # perform permutation
    sig_counts <-
        matrix(0L,
            nrow = length(sc_clust),
            ncol = length(ref_clust)
        )

    for (i in seq_len(n_perm)) {
        resampled <- sample(cluster_ids,
            length(cluster_ids),
            replace = FALSE
        )

        if (!per_cell) {
            permuted_avg <- average_clusters(
                expr_mat,
                resampled,
                method = pseudobulk_method
            )
        } else {
            permuted_avg <- expr_mat[, resampled, drop = FALSE]
        }

        # permutate assignment
        new_score <- calc_similarity(permuted_avg,
            ref_mat,
            compute_method,
            rm0 = rm0,
            ...
        )
        sig_counts <-
            sig_counts + as.numeric(new_score > assigned_score)
    }

    rownames(assigned_score) <- sc_clust
    colnames(assigned_score) <- ref_clust
    rownames(sig_counts) <- sc_clust
    colnames(sig_counts) <- ref_clust

    return(list(
        score = assigned_score,
        p_val = sig_counts / n_perm
    ))
}

#' compute similarity
#' @param query_mat query data matrix
#' @param ref_mat reference data matrix
#' @param compute_method method(s) for computing similarity scores
#' @param rm0 consider 0 as missing data, recommended for per_cell
#' @param ...  additional parameters
#' @return matrix of numeric values
calc_similarity <- function(
    query_mat,
    ref_mat,
    compute_method,
    rm0 = FALSE,
    ...) {
    # remove 0s ?
    if (rm0) {
        message("considering 0 as missing data")
        query_mat[query_mat == 0] <- NA
        similarity_score <- suppressWarnings(stats::cor(as.matrix(query_mat),
            ref_mat,
            method = compute_method,
            use = "pairwise.complete.obs"
        ))
        return(similarity_score)
    } else {
        if (any(compute_method %in% c(
            "pearson",
            "spearman",
            "kendall"
        ))) {
            similarity_score <- suppressWarnings(
                stats::cor(as.matrix(query_mat),
                    ref_mat,
                    method = compute_method
                )
            )
            return(similarity_score)
        }
        if (compute_method == "cosine") {
            res <- proxy::simil(as.matrix(query_mat), ref_mat, method = "cosine", by_rows = FALSE)
            similarity_score <- matrix(res, nrow = nrow(res))
            return(similarity_score)
        }
    }

    sc_clust <- colnames(query_mat)
    ref_clust <- colnames(ref_mat)
    features <- intersect(rownames(query_mat), rownames(ref_mat))
    query_mat <- query_mat[features, ]
    ref_mat <- ref_mat[features, ]
    similarity_score <- matrix(NA,
        nrow = length(sc_clust),
        ncol = length(ref_clust)
    )
    for (i in seq_along(sc_clust)) {
        for (j in seq_along(ref_clust)) {
            similarity_score[i, j] <- vector_similarity(
                query_mat[, sc_clust[i]],
                ref_mat[, ref_clust[j]],
                compute_method, ...
            )
        }
    }
    return(similarity_score)
}

#' Compute similarity between two vectors
#'
#' @description Compute the similarity score between two vectors using a
#' customized scoring function
#' Two vectors may be from either scRNA-seq or bulk RNA-seq data.
#' The lengths of vec1 and vec2 must match, and must be arranged in the
#' same order of genes.
#' Both vectors should be provided to this function after pre-processing,
#' feature selection and dimension reduction.
#'
#' @param vec1 test vector
#' @param vec2 reference vector
#' @param compute_method method to run i.e. corr_coef
#' @param ... arguments to pass to compute_method function
#' @return numeric value of desired correlation or distance measurement
vector_similarity <- function(vec1, vec2, compute_method, ...) {
    # examine whether two vectors are of the same size
    if (!is.numeric(vec1) ||
        !is.numeric(vec2) || length(vec1) != length(vec2)) {
        stop(
            "compute_similarity: two input vectors",
            " are not numeric or of different sizes.",
            call. = FALSE
        )
    }

    if (!(compute_method %in% c("cosine", "kl_divergence"))) {
        stop(compute_method, " not implemented", call. = FALSE)
    }

    if (compute_method == "kl_divergence") {
        res <- kl_divergence(vec1, vec2, ...)
    } else if (compute_method == "cosine") {
        res <- cosine(vec1, vec2, ...)
    }
    # return the similarity score, must be
    return(res)
}

#' Cosine distance
#' @param vec1 test vector
#' @param vec2 reference vector
#' @return numeric value of cosine distance between the vectors
cosine <- function(vec1, vec2) {
    sum(vec1 * vec2) / sqrt(sum(vec1^2) * sum(vec2^2))
}
#' KL divergence
#'
#' @description Use package entropy to compute Kullback-Leibler divergence.
#' The function first converts each vector's reads to pseudo-number of
#' transcripts by normalizing the total reads to total_reads.
#' The normalized read for each gene is then rounded to serve as the
#' pseudo-number of transcripts.
#' Function entropy::KL.shrink is called to compute the KL-divergence between
#' the two vectors, and the maximal allowed divergence is set to max_KL.
#' Finally, a linear transform is performed to convert the KL divergence,
#' which is between 0 and max_KL, to a similarity score between -1 and 1.
#'
#' @param vec1 Test vector
#' @param vec2 Reference vector
#' @param if_log Whether the vectors are log-transformed. If so, the
#' raw count should be computed before computing KL-divergence.
#' @param total_reads Pseudo-library size
#' @param max_KL Maximal allowed value of KL-divergence.
#' @return numeric value, with additional attributes, of kl divergence
#'  between the vectors
kl_divergence <- function(
    vec1,
    vec2,
    if_log = FALSE,
    total_reads = 1000,
    max_KL = 1) {
    if (if_log) {
        vec1 <- expm1(vec1)
        vec2 <- expm1(vec2)
    }
    count1 <- round(vec1 * total_reads / sum(vec1))
    count2 <- round(vec2 * total_reads / sum(vec2))
    est_KL <- entropy::KL.shrink(count1, count2,
        unit = "log",
        verbose = FALSE
    )
    return((max_KL - est_KL) / max_KL * 2 - 1)
}
