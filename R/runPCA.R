#' Perform Principal components analysis with `scran.chan`
#'
#' @param object A matrix-like object containing the counts (non-negative
#' integers). Rows are features and columns are cells.
#' @param ... Additional arguments passed to default methods.
#' @export
runPCA <- function(object, ...) UseMethod("runPCA")

#' @param exprs_values Integer scalar or string indicating which assay of x
#' contains the expression values.
#' @param dimred String or integer scalar specifying the existing dimensionality
#' reduction results to use.
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if
#' dimred is specified.
#' @param size_factors A numeric vector of length equal to the number of cells
#' in x, containing positive size factors for all cells.
#' @param name String specifying the name to be used to store the result in the
#' [reducedDims][SingleCellExperiment::reducedDim] of the output.
#' @export
#' @rdname runPCA
runPCA.SingleCellExperiment <- function(object, ...,
                                        exprs_values = "counts",
                                        dimred = NULL, n_dimred = NULL,
                                        size_factors = NULL,
                                        name = "PCA") {
    size_factors <- size_factors %||% SingleCellExperiment::sizeFactors(object)
    mat <- .get_mat_from_sce(object, exprs_values, dimred, n_dimred)
    pca <- runPCA(object = mat, ..., size_factors = size_factors)
    SingleCellExperiment::reducedDim(object, name) <- pca
    object
}

#' @param d Integer scalar specifying the number of top PCs to obtain.
#' @param scale Logical scalar indicating whether to scale rows to unit
#' variance.
#' @param subset_row Integer, logical or character vector specifying which
#' features to use in the PCA (e.g., highly variable genes). If NULL, all
#' features in x are used.
#' @param batch Vector or factor of length equal to the number of cells,
#' specifying the batch of origin for each cell. Alternatively NULL if all cells
#' belong to the same batch.
#' @param norm_batch String indicating how batch should be handled when
#' centering the size factors. If `"lowest"`, we downscale all batches to the
#' coverage of the lowest batch. If `"perblock"`, we scale each batch to a mean
#' of 1.
#' @param pca_batch String indicating how batch should be handled (if it is
#' supplied). "block" is equivalent to linear regression on x prior to PCA,
#' while "weight" will only weight each batch so that they contribute equally to
#' the PCA.
#' @param force_integer Logical scalar indicating whether double-precision x
#' should be forced into integers.
#' @param no_sparse_copy Logical scalar indicating whether we should avoid a
#' copy when object is a [dgCMatrix][Matrix::dgCMatrix-class] This is more
#' memory efficient if the data has already been loaded into memory. If TRUE,
#' any setting of force.integer is ignored.
#' @param threads Integer scalar specifying the number of threads to use.
#' @seealso [runPCA.chan][scran.chan::runPCA.chan]
#' @export
#' @rdname runPCA
runPCA.default <- function(object, d = 50L, scale = FALSE, ...,
                           size_factors = NULL, subset_row = NULL,
                           batch = NULL, norm_batch = NULL, pca_batch = NULL,
                           force_integer = TRUE, no_sparse_copy = TRUE,
                           threads = 1L) {
    # nromalization, adjust for differences in sequencing depth --------
    norm <- scran.chan::logNormCounts.chan(
        x = scran.chan::initializeSparseMatrix(
            object,
            force.integer = force_integer,
            no.sparse.copy = no_sparse_copy,
            by.column = TRUE,
            num.threads = threads
        ),
        size.factors = size_factors,
        batch = batch,
        batch.mode = norm_batch,
        num.threads = threads
    )

    # dimensionality reduction -----------------------------------------
    batch_pcs <- scran.chan::runPCA.chan(
        x = norm,
        num.comp = as.integer(d),
        subset = subset_row,
        scale = scale,
        num.threads = threads,
        batch = batch,
        batch.method = pca_batch,
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
