#' Fast mutual nearest neighbors correction with `scran.chan`
#'
#' @export
runMNN <- function(object, ...) UseMethod("runMNN")

#' @inheritParams runPCA
#' @export
#' @rdname runMNN
runMNN.SingleCellExperiment <- function(object, ...,
                                        dimred = "PCA", n_dimred = NULL,
                                        exprs_values = NULL,
                                        name = "corrected") {
    mat <- .get_mat_from_sce(object, exprs_values, dimred, n_dimred)
    mnn <- runMNN(object = mat, ...)
    SingleCellExperiment::reducedDim(object, name) <- mnn
    object
}

#' @param k Integer scalar specifying the number of neighbors to use when
#' identifying MNN pairs.
#' @inheritDotParams scran.chan::mnnCorrect.chan -x -batch -k -approximate -num.threads
#' @param approximate Logical scalar specifying whether to perform an
#' approximate neighbor search.
#' @seealso [mnnCorrect.chan][scran.chan::mnnCorrect.chan]
#' @export
#' @rdname runMNN
runMNN.default <- function(object, batch = NULL, k = 15L,
                           ..., approximate = TRUE, threads = 1L) {
    assert_bool(approximate)
    threads <- as.integer(threads)
    # run MNN --------------------------------------------------------
    mnn_res <- scran.chan::mnnCorrect.chan(
        x = object,
        batch = batch, k = as.integer(k), ...,
        approximate = approximate,
        num.threads = threads
    )
    out <- t(mnn_res$corrected)
    attr(out, "merge.order") <- mnn_res$merge.order
    attr(out, "num.pairs") <- mnn_res$num.pairs
    out
}
