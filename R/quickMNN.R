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
    batch, ...,
    # MNN arguments
    k = 15L, nmads = 3L, mass_cap = NULL,
    order = NULL, reference_policy = NULL,
    approximate = TRUE, threads = NULL, name = "corrected") {
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
        approximate = approximate, threads = threads,
        name = name
    )
}

#' @inheritParams runMNN
#' @inheritDotParams runPCA
#' @export
#' @rdname quickMNN
quickMNN.SingleCellExperiment <- .quickMNN

#' @export
#' @rdname quickMNN
quickMNN.Seurat <- .quickMNN
