#' Score feature set activity for each cell
#'
#' @inheritParams runPCA
#' @export
scoreFeatureSet <- function(object, ...) UseMethod("scoreFeatureSet")

#' @export
#' @rdname scoreFeatureSet
scoreFeatureSet.SingleCellExperiment <- function(object, ...,
                                                 assay = "logcounts") {
    mat <- .get_mat_from_sce(object, assay)
    scoreFeatureSet(object = mat, ...)
}

#' @export
#' @rdname scoreFeatureSet
scoreFeatureSet.Seurat <- function(object, ..., assay = NULL, layer = "data") {
    mat <- .get_mat_from_seurat(object, assay, layer)
    scoreFeatureSet(object = mat, ...)
}

#' @inheritParams scran.chan::scoreFeatureSet.chan
#' @inherit scran.chan::scoreFeatureSet.chan description details return
#' @seealso [scoreFeatureSet.chan][scran.chan::scoreFeatureSet.chan]
#' @rdname scoreFeatureSet
#' @export
scoreFeatureSet.default <- function(object, features,
                                    batch = NULL, scale = FALSE,
                                    ...,
                                    force_integer = FALSE,
                                    no_sparse_copy = TRUE,
                                    threads = NULL) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    scran.chan::scoreFeatureSet.chan(
        x = initialize_matrix(
            object,
            force_integer = force_integer,
            no_sparse_copy = no_sparse_copy,
            threads = threads
        ),
        features = features, batch = batch, scale = scale,
        num.threads = threads
    )
}
