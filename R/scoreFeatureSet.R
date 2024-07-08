#' Score feature set activity for each cell
#'
#' @inheritParams runPCA
#' @export
scoreFeatureSet <- function(object, ...) UseMethod("scoreFeatureSet")

#' @param name A string of score assay name.
#' @export
#' @rdname scoreFeatureSet
scoreFeatureSet.SummarizedExperiment <- function(object, ...,
                                                 name = "scores",
                                                 assay = "logcounts") {
    mat <- .get_mat_from_sce(object, assay)
    scores <- list(scoreFeatureSet(object = mat, ...))
    names(scores) <- name
    SummarizedExperiment::SummarizedExperiment(
        assays = scores,
        colData = SummarizedExperiment::colData(object)
    )
}

#' @export
#' @rdname scoreFeatureSet
scoreFeatureSet.Seurat <- function(object, ..., assay = NULL, layer = "data") {
    mat <- .get_mat_from_seurat(object, assay, layer)
    scoreFeatureSet(object = mat, ...)
}

#' @param feature_sets A list of Integer, logical or character vector specifying
#' the features that belong to the set. You can also directly input an atomic
#' vector if there is only one feature set.
#' @inheritParams scran.chan::scoreFeatureSet.chan
#' @inherit scran.chan::scoreFeatureSet.chan description details
#' @return A matrix of feature scores.
#' @seealso [scoreFeatureSet.chan][scran.chan::scoreFeatureSet.chan]
#' @rdname scoreFeatureSet
#' @export
scoreFeatureSet.default <- function(object, feature_sets,
                                    batch = NULL, scale = FALSE,
                                    ...,
                                    force_integer = FALSE,
                                    no_sparse_copy = TRUE,
                                    threads = NULL) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    mat <- initialize_matrix(
        object,
        force_integer = force_integer,
        no_sparse_copy = no_sparse_copy,
        threads = threads
    )
    if (is.list(feature_sets)) {
        out <- lapply(feature_sets, function(features) {
            .scoreFeatureSet(mat = mat, features, batch, scale, threads)
        })
        weights <- lapply(out, attr, which = "weights", exact = TRUE)
        names(weights) <- paste(
            names(weights) %||% seq_along(weights),
            "weights",
            sep = "_"
        )
        out <- do.call(base::rbind, out)
        colnames(out) <- colnames(object)
        rlang::inject(structure(out, !!!weights))
    } else {
        out <- .scoreFeatureSet(mat = mat, feature_sets, batch, scale, threads)
        dim(out) <- c(1L, length(out))
        colnames(out) <- colnames(object)
        out
    }
}

.scoreFeatureSet <- function(mat, features, batch, scale, threads) {
    out <- scran.chan::scoreFeatureSet.chan(
        x = mat,
        features = features, batch = batch, scale = scale,
        num.threads = threads
    )
    structure(out$scores, weights = out$weights)
}
