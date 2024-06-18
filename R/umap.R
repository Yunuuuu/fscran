#' Fast mutual nearest neighbors correction with `scran.chan`
#'
#' @seealso [mnnCorrect.chan][scran.chan::mnnCorrect.chan]
#' @export
chan_umap <- function(object, ...) UseMethod("chan_umap")

#' @inheritParams chan_umap
#' @export
#' @rdname chan_mnn
chan_umap.SingleCellExperiment <- function(object, ...,
                                           exprs_values = "counts",
                                           dimred = NULL, n_dimred = NULL,
                                           name = "UMAP") {
    mat <- .get_mat_from_sce(object, exprs_values, dimred, n_dimred)
    umap <- chan_umap(object = mat, ...)
    SingleCellExperiment::reducedDim(object, name) <- umap
    object
}

#' @param n_neighbors Integer scalar specifying the number of neighbors to use
#' in the UMAP algorithm.
#' @param min_dist Numeric scalar specifying the minimum distance between
#' points.
#' @param n_epochs Integer scalar specifying the number of epochs to perform. If
#' set to `-1`, an appropriate number of epochs is chosen based on ncol(x).
#' @param seed Integer scalar specifying the seed to use.
#' @inheritParams chan_pca
#' @inheritParams scran.chan::runUMAP.chan
#' @export
#' @rdname chan_umap
chan_umap.default <- function(object, n_neighbors = 20L,
                              min_dist = 0.01, n_epochs = -1L,
                              ...,
                              approximate = TRUE, seed = 123456L,
                              threads = 1L) {
    assert_number(n_neighbors)
    assert_number(min_dist)
    assert_number(n_epochs)
    assert_number(seed)
    scran.chan::runUMAP.chan(
        x = object,
        num.neighbors = as.integer(n_neighbors),
        num.epochs = n_epochs,
        min.dist = min_dist,
        seed = as.integer(seed),
        approximate = approximate,
        downsample = NULL,
        num.threads = as.integer(threads)
    )
}
