#' Perform Principal components analysis with `scran.chan`
#'
#' @inheritParams runPCA
#' @inherit runPCA return
#' @seealso
#' - [logNormCounts]
#' - [runPCA]
#' @export
logNormAndPCA <- function(object, ...) UseMethod("logNormAndPCA")

#' @export
#' @rdname logNormAndPCA
logNormAndPCA.SingleCellExperiment <- function(object, ..., name = "PCA") {
    pca <- NextMethod("logNormAndPCA", object = object, ...)
    add_dimred_to_sce(object, pca, name)
}

#' @export
#' @rdname logNormAndPCA
logNormAndPCA.Seurat <- function(object, ...,
                                 assay = NULL, layer = "counts",
                                 name = "PCA") {
    pca <- NextMethod("logNormAndPCA",
        object = object, ..., assay = assay, layer = layer
    )
    add_dimred_to_seurat(object, pca, name, assay, layer, NULL)
}

#' @inheritDotParams logNormCounts
#' @export
#' @rdname logNormAndPCA
logNormAndPCA.default <- function(object, d = 50L, scale = FALSE, ...,
                                  subset_row = NULL,
                                  batch = NULL, batch_method = NULL) {
    threads <- set_threads(threads)
    # nromalization, adjust for differences in sequencing depth --------
    logcounts <- logNormCounts(
        object = object,
        batch = batch,
        num.threads = threads,
        ...
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
