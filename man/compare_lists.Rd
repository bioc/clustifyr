% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_genelist.R
\name{compare_lists}
\alias{compare_lists}
\title{Calculate adjusted p-values for hypergeometric test of gene lists
or jaccard index}
\usage{
compare_lists(
  bin_mat,
  marker_mat,
  n = 30000,
  metric = "hyper",
  output_high = TRUE,
  details_out = FALSE
)
}
\arguments{
\item{bin_mat}{binarized single-cell expression matrix,
feed in by_cluster mat, if desired}

\item{marker_mat}{matrix or dataframe of candidate genes for each cluster}

\item{n}{number of genes in the genome}

\item{metric}{adjusted p-value for hypergeometric test, or jaccard index}

\item{output_high}{if true (by default to fit with rest of package),
-log10 transform p-value}

\item{details_out}{whether to also output shared gene list from jaccard}
}
\value{
matrix of numeric values, clusters from expr_mat as row names,
cell types from marker_mat as column names
}
\description{
Calculate adjusted p-values for hypergeometric test of gene lists
or jaccard index
}
\examples{
pbmc_mm <- matrixize_markers(pbmc_markers)

pbmc_avg <- average_clusters(
    pbmc_matrix_small,
    pbmc_meta,
    cluster_col = "classified"
)

pbmc_avgb <- binarize_expr(pbmc_avg)

compare_lists(
    pbmc_avgb,
    pbmc_mm,
    metric = "spearman"
)
}
