#' Compute the uniform manifold approximation and projection with `scran.chan`
#'
#' @export
runUMAP <- function(object, ...) UseMethod("runUMAP")

#' @inheritParams runPCA
#' @export
#' @rdname runUMAP
runUMAP.SingleCellExperiment <- function(object, ...,
                                         dimred = "PCA", n_dimred = NULL,
                                         exprs_values = NULL,
                                         name = "UMAP") {
    mat <- .get_mat_from_sce(object, exprs_values, dimred, n_dimred)
    umap <- runUMAP(object = mat, ...)
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
#' @inheritParams runPCA
#' @inheritParams scran.chan::runUMAP.chan
#' @inherit runPCA return
#' @seealso [runUMAP.chan][scran.chan::runUMAP.chan]
#' @export
#' @rdname runUMAP
runUMAP.default <- function(object, n_neighbors = 15L,
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
