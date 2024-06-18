#' Perform Principal components analysis with `scran.chan`
#'
#' @param object A matrix-like object. Rows are features and columns are cells.
#' @inheritParams logNormCounts
#' @export
runPCA <- function(object, ...) UseMethod("runPCA")

#' @param dimred String or integer scalar specifying the existing dimensionality
#' reduction results to use.
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if
#' dimred is specified.
#' @param name String specifying the name to be used to store the result in the
#' [reducedDims][SingleCellExperiment::reducedDim] or
#' [reductions][SeuratObject::Seurat-class] of the output.
#' @export
#' @rdname runPCA
runPCA.SingleCellExperiment <- function(object, ...,
                                        assay = "logcounts",
                                        dimred = NULL, n_dimred = NULL,
                                        name = "PCA") {
    mat <- .get_mat_from_sce(object, assay, dimred, n_dimred)
    pca <- runPCA(object = mat, ...)
    add_dimred_to_sce(object, pca, name)
}

#' @export
#' @rdname runPCA
runPCA.Seurat <- function(object, ...,
                          assay = NULL, layer = "data",
                          dimred = NULL, n_dimred = NULL,
                          name = "PCA") {
    mat <- .get_mat_from_seurat(object, assay, layer, dimred, n_dimred)
    pca <- runPCA(object = mat, ...)
    add_dimred_to_seurat(object, pca, name, assay, layer, dimred)
}

#' @param d Integer scalar specifying the number of top PCs to obtain.
#' @param scale Logical scalar indicating whether to scale rows to unit
#' variance.
#' @param subset_row Integer, logical or character vector specifying which
#' features to use in the PCA (e.g., highly variable genes). If NULL, all
#' features in x are used.
#' @param batch_method String indicating how batch should be handled (if it is
#' supplied). `"block"` is equivalent to linear regression on x prior to PCA,
#' while `"weight"` will only weight each batch so that they contribute equally
#' to the PCA.
#' @inheritParams logNormCounts
#' @seealso [runPCA.chan][scran.chan::runPCA.chan]
#' @return
#'  - `default` method: A numeric matrix where rows are cells and columns are
#'                    the two dimensions of the embedding.
#'  - `SingleCellExperiment` method: embedding was added into
#'    [reducedDims][SingleCellExperiment::reducedDim] named as `name`.
#'  - `Seurat` method: embedding was added into
#'    [reductions][SeuratObject::Seurat-class] named as `name`.
#' @export
#' @rdname runPCA
runPCA.default <- function(object, d = 50L, scale = FALSE, ...,
                           subset_row = NULL, batch = NULL, batch_method = NULL,
                           force_integer = TRUE, no_sparse_copy = TRUE,
                           threads = NULL) {
    .runPCA(
        logcounts = scran.chan::initializeSparseMatrix(
            object,
            force.integer = force_integer,
            no.sparse.copy = no_sparse_copy,
            by.column = TRUE,
            num.threads = threads
        ),
        d = d,
        subset_row = subset_row,
        scale = scale,
        batch = batch,
        batch_method = batch_method,
        threads = set_threads(threads)
    )
}

.runPCA <- function(logcounts, d, subset_row, scale,
                    batch, batch_method, threads) {
    batch_pcs <- scran.chan::runPCA.chan(
        x = logcounts,
        num.comp = as.integer(d),
        subset = subset_row,
        scale = scale,
        num.threads = threads,
        batch = batch,
        batch.method = batch_method,
        rotation = TRUE
    )
    # batch_pcs$components:
    # a numeric matrix containing the top principal components. Each row
    # corresponds to a PC and each column corresponds to a cell.
    out <- t(batch_pcs$components)
    attr(out, "percentVar") <- batch_pcs$prop.variance
    attr(out, "rotation") <- batch_pcs$rotation
    out
}
