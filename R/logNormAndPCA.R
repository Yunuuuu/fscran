#' Perform Principal components analysis with `scran.chan`
#'
#' @inheritParams logNormCounts
#' @inheritParams runPCA
#' @inherit runPCA return
#' @seealso 
#' - [logNormCounts]
#' - [runPCA]
#' @export
logNormAndPCA <- function(object, ...) UseMethod("logNormAndPCA")


#' @export
#' @rdname logNormAndPCA
logNormAndPCA.SingleCellExperiment <- function(object, ...,
                                               size_factors = NULL,
                                               assay = "counts",
                                               name = "PCA") {
    size_factors <- size_factors %||% SingleCellExperiment::sizeFactors(object)
    mat <- .get_mat_from_sce(object, assay, NULL, NULL)
    pca <- logNormAndPCA(object = mat, ..., size_factors = size_factors)
    add_dimred_to_sce(object, pca, name)
}

#' @export
#' @rdname logNormAndPCA
logNormAndPCA.Seurat <- function(object, ...,
                                 assay = NULL, layer = "counts",
                                 name = "PCA") {
    mat <- .get_mat_from_seurat(object, assay, layer, NULL, NULL)
    pca <- logNormAndPCA(object = mat, ...)
    add_dimred_to_seurat(object, pca, name, assay, layer, NULL)
}

#' @export
#' @rdname logNormAndPCA
logNormAndPCA.default <- function(object, d = 50L, scale = FALSE, ...,
                                  size_factors = NULL, subset_row = NULL,
                                  batch = NULL, batch_mode = NULL,
                                  batch_method = NULL,
                                  force_integer = TRUE, no_sparse_copy = TRUE,
                                  threads = NULL) {
    batch_mode <- match.arg(batch_mode, c("perblock", "lowest"))
    threads <- set_threads(threads)
    # nromalization, adjust for differences in sequencing depth --------
    logcounts <- scran.chan::logNormCounts.chan(
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

    # dimensionality reduction -----------------------------------------
    .runPCA(
        logcounts = logcounts,
        d = d,
        subset_row = subset_row,
        scale = scale,
        batch = batch,
        batch_method = batch_method,
        threads = threads
    )
}
