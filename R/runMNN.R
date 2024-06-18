#' Fast mutual nearest neighbors correction with `scran.chan`
#'
#' @export
runMNN <- function(object, ...) UseMethod("runMNN")

#' @inheritParams runPCA
#' @export
#' @rdname runMNN
runMNN.SingleCellExperiment <- function(object, ...,
                                        dimred = "PCA", n_dimred = NULL,
                                        assay = NULL,
                                        name = "corrected") {
    mat <- .get_mat_from_sce(object, assay, dimred, n_dimred)
    mnn <- runMNN(object = mat, ...)
    add_dimred_to_sce(object, mnn, name)
}

#' @inheritParams runPCA
#' @export
#' @rdname runMNN
runMNN.Seurat <- function(object, ...,
                          dimred = "PCA", n_dimred = NULL,
                          assay = NULL, layer = NULL,
                          name = "corrected") {
    mat <- .get_mat_from_seurat(object, assay, layer, dimred, n_dimred)
    mnn <- runMNN(object = mat, ...)
    add_dimred_to_seurat(object, mnn, name, assay, layer, dimred)
}

#' @param k Integer scalar specifying the number of neighbors to use when
#' identifying MNN pairs.
#' @inheritParams scran.chan::mnnCorrect.chan
#' @param mass_cap Integer scalar specifying the cap on the number of
#' observations to use for center of mass calculations. A value of 100,000 may
#' be appropriate for speeding up correction of very large datasets. If this is
#' set to `NULL`, no cap is used.
#' @param reference_policy String specifying the policy to use to choose the
#' first reference batch. This can be based on the largest batch (`"max-size"`),
#' the most variable batch (`"max-variance"`), some combination of those two
#' (`"max-rss"`) or the first specified input ("input"). Only used for automatic
#' merges, i.e., when order=`NULL`.
#' @param approximate Logical scalar specifying whether to perform an
#' approximate neighbor search.
#' @seealso [mnnCorrect.chan][scran.chan::mnnCorrect.chan]
#' @export
#' @rdname runMNN
runMNN.default <- function(object, batch, k = 15L, ...,
                           nmads = 3L, mass_cap = NULL,
                           order = NULL, reference_policy = NULL,
                           approximate = TRUE, threads = NULL) {
    # run MNN --------------------------------------------------------
    .runMNN(
        object = object,
        batch = batch, k = k,
        nmads = nmads, mass_cap = mass_cap,
        order = order,
        reference_policy = reference_policy,
        approximate = approximate,
        threads = set_threads(threads)
    )
}

.runMNN <- function(object, batch, k, nmads, mass_cap,
                    order, reference_policy, approximate,
                    threads) {
    mnn_res <- scran.chan::mnnCorrect.chan(
        x = object,
        batch = batch, k = as.integer(k),
        nmads = nmads, mass.cap = mass_cap,
        order = order, reference.policy = reference_policy,
        approximate = approximate,
        num.threads = threads
    )
    out <- t(mnn_res$corrected)
    attr(out, "merge.order") <- mnn_res$merge.order
    attr(out, "num.pairs") <- mnn_res$num.pairs
    out
}
