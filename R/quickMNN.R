#' Quick fastMNN with `scran.chan`
#'
#' @seealso
#' - [logNormCounts.chan][scran.chan::logNormCounts.chan]
#' - [runPCA.chan][scran.chan::runPCA.chan]
#' - [mnnCorrect.chan][scran.chan::mnnCorrect.chan]
#' @export
quickMNN <- function(object, ...) UseMethod("quickMNN")

.quickMNN <- function(
    object, # PCA arguments
    batch = NULL, ...,
    # MNN arguments
    k = 15L, nmads = 3L, mass_cap = NULL,
    order = NULL, reference_policy = NULL,
    approximate = TRUE, threads = 1L) {
    # dimensionality reduction --------------
    object <- runPCA(
        object = object, batch = batch, ...,
        threads = threads, name = "PCA"
    )
    # run MNN -------------------------------
    runMNN(
        object = object, batch = batch, dimred = "PCA",
        k = k, nmads = nmads, mass_cap = mass_cap,
        order = order, reference_policy = reference_policy,
        approximate = approximate, threads = threads
    )
}

#' @inheritParams runMNN
#' @inheritDotParams runPCA
#' @export
#' @rdname quickMNN
quickMNN.SingleCellExperiment <- .quickMNN
