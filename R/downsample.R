#' Downsample cells based on their neighbors
#'
#' @export
downsample <- function(object, ...) UseMethod("downsample")

#' @inheritParams runPCA
#' @export
#' @rdname downsample
downsample.SingleCellExperiment <- function(object, ...,
                                            dimred = "PCA", n_dimred = NULL,
                                            assay = NULL) {
    mat <- .get_mat_from_sce(object, assay, dimred, n_dimred)
    downsample(object = mat, ...)
}

#' @export
#' @rdname downsample
downsample.SingleCellExperiment <- function(object, ...,
                                            dimred = "PCA", n_dimred = NULL,
                                            assay = NULL, layer = NULL) {
    mat <- .get_mat_from_seurat(object, assay, layer, dimred, n_dimred)
    downsample(object = mat, ...)
}

#' @inheritParams runPCA
#' @inheritDotParams scran.chan::downsampleByNeighbors.chan -x -num.threads
#' @seealso [downsampleByNeighbors.chan][scran.chan::downsampleByNeighbors.chan]
#' @return List containing:
#' - `chosen`: an integer vector containing the column indices of x to retain.
#' - `assigned`: integer vector of length equal to `ncol(object)`, containing
#'   the column index of the representative cell for each cell in `object`. 
#' @export
#' @rdname downsample
downsample.default <- function(object, ..., threads = NULL) {
    scran.chan::downsampleByNeighbors.chan(
        x = object, ...,
        num.threads = set_threads(threads)
    )
}
