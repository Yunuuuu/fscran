#' Log-normalize the counts
#'
#' @inheritParams runPCA
#' @export
logNormCounts <- function(object, ...) UseMethod("logNormCounts")

#' @export
#' @rdname logNormCounts
logNormCounts.SingleCellExperiment <- function(object, ...,
                                               size_factors = NULL,
                                               assay = "counts") {
    size_factors <- size_factors %||% SingleCellExperiment::sizeFactors(object)
    mat <- .get_mat_from_sce(object, assay, NULL, NULL)
    logNormCounts(object = mat, size_factors = size_factors, ...)
}

#' @export
#' @rdname logNormCounts
logNormCounts.Seurat <- function(object, ...,
                                 assay = NULL, layer = "counts") {
    mat <- .get_mat_from_seurat(object, assay, layer, NULL, NULL)
    logNormCounts(object = mat, ...)
}

#' @param size_factors A numeric vector of length equal to the number of cells
#' in x, containing positive size factors for all cells.
#' @param batch_mode String indicating how batch should be handled when
#' centering the size factors. If `"lowest"`, we downscale all batches to the
#' coverage of the lowest batch. If `"perblock"`, we scale each batch to a mean
#' of 1. Default: `"perblock"`.
#' @inherit scran.chan::logNormCounts.chan return
#' @export
#' @rdname logNormCounts
logNormCounts.default <- function(object, ...,
                                  size_factors = NULL,
                                  batch = NULL, batch_mode = NULL,
                                  force_integer = TRUE, no_sparse_copy = TRUE,
                                  threads = NULL) {
    batch_mode <- match.arg(batch_mode, c("perblock", "lowest"))
    threads <- set_threads(threads)
    scran.chan::logNormCounts.chan(
        x = scran.chan::initializeSparseMatrix(
            object,
            force.integer = force_integer,
            no.sparse.copy = no_sparse_copy,
            by.column = TRUE,
            num.threads = threads
        ),
        size.factors = size_factors,
        batch = batch,
        batch.mode = batch_mode,
        num.threads = threads
    )
}
