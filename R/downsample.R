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

#' @inheritParams runPCA
#' @inheritDotParams scran.chan::downsampleByNeighbors.chan -x -num.threads
#' @seealso [downsampleByNeighbors.chan][scran.chan::downsampleByNeighbors.chan]
#' @inherit scran.chan::downsampleByNeighbors.chan return
#' @export
#' @rdname downsample
downsample.default <- function(object, ..., threads = 1L) {
    scran.chan::downsampleByNeighbors.chan(
        x = object, ...,
        num.threads = threads
    )
}
