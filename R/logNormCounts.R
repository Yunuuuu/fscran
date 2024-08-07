#' Log-normalize the counts
#'
#' @param object A matrix-like object containing the counts (non-negative
#' integers). Rows are features and columns are cells.
#' @param ...
#'  - `default` method: cannot be used.
#'  - `SingleCellExperiment` and `Seurat` method: additional arguments passed to
#'    `default` methods.
#' @export
logNormCounts <- function(object, ...) UseMethod("logNormCounts")

#' @param assay Integer scalar or string indicating which assay of x
#' contains the expression values.
#' @export
#' @rdname logNormCounts
logNormCounts.SingleCellExperiment <- function(object, size_factors = NULL,
                                               ..., assay = "counts") {
    size_factors <- size_factors %||% SingleCellExperiment::sizeFactors(object)
    mat <- .get_mat_from_sce(object, assay)
    logNormCounts(object = mat, size_factors = size_factors, ...)
}

#' @param layer Name of the layer to get from the assay data.
#' @export
#' @rdname logNormCounts
logNormCounts.Seurat <- function(object, ..., assay = NULL, layer = "counts") {
    mat <- .get_mat_from_seurat(object, assay, layer)
    logNormCounts(object = mat, ...)
}

#' @param size_factors A numeric vector of length equal to the number of cells
#' in x, containing positive size factors for all cells.
#' @param batch Vector or factor of length equal to the number of cells,
#' specifying the batch of origin for each cell. Alternatively `NULL` if all
#' cells belong to the same batch.
#' @param batch_mode String indicating how batch should be handled when
#' centering the size factors. If `"lowest"`, we downscale all batches to the
#' coverage of the lowest batch. If `"perblock"`, we scale each batch to a mean
#' of 1. Default: `"perblock"`.
#' @param force_integer Logical scalar indicating whether double-precision
#' `object` should be forced into integers.
#' @param no_sparse_copy Logical scalar indicating whether we should avoid a
#' copy when object is a [dgCMatrix][Matrix::dgCMatrix-class] This is more
#' memory efficient if the data has already been loaded into memory. If `TRUE`,
#' any setting of force.integer is ignored.
#' @param threads Integer scalar specifying the number of threads to use. If
#' `NULL`, all detected threads will be used. See
#' [detectCores][parallel::detectCores].
#' @inherit scran.chan::logNormCounts.chan return
#' @seealso [logNormCounts.chan][scran.chan::logNormCounts.chan]
#' @export
#' @rdname logNormCounts
logNormCounts.default <- function(object, size_factors = NULL,
                                  batch = NULL, batch_mode = NULL, ...,
                                  force_integer = TRUE, no_sparse_copy = TRUE,
                                  threads = NULL) {
    rlang::check_dots_empty()
    batch_mode <- match.arg(batch_mode, c("perblock", "lowest"))
    threads <- set_threads(threads)
    scran.chan::logNormCounts.chan(
        x = initialize_matrix(
            object,
            force_integer = force_integer,
            no_sparse_copy = no_sparse_copy,
            threads = threads
        ),
        size.factors = size_factors,
        batch = batch,
        batch.mode = batch_mode,
        num.threads = threads
    )
}
