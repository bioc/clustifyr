% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{clustify_nudge}
\alias{clustify_nudge}
\alias{clustify_nudge.default}
\alias{clustify_nudge.Seurat}
\title{Combined function to compare scRNA-seq data to
bulk RNA-seq data and marker list}
\usage{
clustify_nudge(input, ...)

\method{clustify_nudge}{default}(
  input,
  ref_mat,
  marker,
  metadata = NULL,
  cluster_col = NULL,
  query_genes = NULL,
  compute_method = "spearman",
  weight = 1,
  threshold = -Inf,
  dr = "umap",
  norm = "diff",
  call = TRUE,
  marker_inmatrix = TRUE,
  mode = "rank",
  obj_out = FALSE,
  seurat_out = obj_out,
  rename_prefix = NULL,
  lookuptable = NULL,
  ...
)

\method{clustify_nudge}{Seurat}(
  input,
  ref_mat,
  marker,
  cluster_col = NULL,
  query_genes = NULL,
  compute_method = "spearman",
  weight = 1,
  obj_out = TRUE,
  seurat_out = obj_out,
  threshold = -Inf,
  dr = "umap",
  norm = "diff",
  marker_inmatrix = TRUE,
  mode = "rank",
  rename_prefix = NULL,
  ...
)
}
\arguments{
\item{input}{express matrix or object}

\item{...}{passed to matrixize_markers}

\item{ref_mat}{reference expression matrix}

\item{marker}{matrix of markers}

\item{metadata}{cell cluster assignments, supplied as a vector
or data.frame. If
data.frame is supplied then \code{cluster_col} needs to be set.}

\item{cluster_col}{column in metadata that contains cluster ids per cell.
Will default to first
column of metadata if not supplied.
Not required if running correlation per cell.}

\item{query_genes}{A vector of genes of interest to compare.
If NULL, then common genes between
the expr_mat and ref_mat will be used for comparision.}

\item{compute_method}{method(s) for computing similarity scores}

\item{weight}{relative weight for the gene list scores,
when added to correlation score}

\item{threshold}{identity calling minimum score threshold,
only used when obj_out = T}

\item{dr}{stored dimension reduction}

\item{norm}{whether and how the results are normalized}

\item{call}{make call or just return score matrix}

\item{marker_inmatrix}{whether markers genes are already
in preprocessed matrix form}

\item{mode}{use marker expression pct or ranked cor score for nudging}

\item{obj_out}{whether to output object instead of cor matrix}

\item{seurat_out}{output cor matrix or called seurat object (deprecated, use obj_out)}

\item{rename_prefix}{prefix to add to type and r column names}

\item{lookuptable}{if not supplied, will look in built-in
table for object parsing}
}
\value{
single cell object, or matrix of numeric values,
clusters from input as row names, cell types from ref_mat as column names
}
\description{
Combined function to compare scRNA-seq data to
bulk RNA-seq data and marker list
}
\examples{

# Seurat
so <- so_pbmc()
clustify_nudge(
    input = so,
    ref_mat = cbmc_ref,
    marker = cbmc_m,
    cluster_col = "seurat_clusters",
    threshold = 0.8,
    obj_out = FALSE,
    mode = "pct",
    dr = "umap"
)

# Matrix
clustify_nudge(
    input = pbmc_matrix_small,
    ref_mat = cbmc_ref,
    metadata = pbmc_meta,
    marker = as.matrix(cbmc_m),
    query_genes = pbmc_vargenes,
    cluster_col = "classified",
    threshold = 0.8,
    call = FALSE,
    marker_inmatrix = FALSE,
    mode = "pct"
)
}
