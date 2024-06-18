#' Compute the t-stochastic neighbor embedding with `scran.chan`
#'
#' @export
runTSNE <- function(object, ...) UseMethod("runTSNE")

#' @inheritParams runPCA
#' @export
#' @rdname runTSNE
runTSNE.SingleCellExperiment <- function(object, ...,
                                         dimred = "PCA", n_dimred = NULL,
                                         assay = NULL,
                                         name = "TSNE") {
    mat <- .get_mat_from_sce(object, assay, dimred, n_dimred)
    tsne <- runTSNE(object = mat, ...)
    SingleCellExperiment::reducedDim(object, name) <- tsne
    object
}

#' @inheritParams runPCA
#' @export
#' @rdname runTSNE
runTSNE.Seurat <- function(object, ...,
                           assay = NULL, layer = "counts",
                           dimred = NULL, n_dimred = NULL,
                           name = "TSNE") {
    mat <- .get_mat_from_seurat(object, assay, layer, dimred, n_dimred)
    umap <- runTSNE(object = mat, ...)
    reduction_key <- SeuratObject::Key(name, quiet = TRUE)
    rownames(umap) <- rownames(mat)
    colnames(umap) <- paste0(reduction_key, seq_len(ncol(umap)))
    object[[name]] <- SeuratObject::CreateDimReducObject(
        embeddings = umap,
        stdev = as.numeric(apply(umap, 2L, stats::sd)),
        assay = .get_assay_from_seurat(object, assay, layer, dimred),
        key = reduction_key
    )
    object
}

#' @param perplexity Numeric scalar specifying the perplexity to use in the
#' t-SNE algorithm.
#' @param interpolate Integer scalar specifying the grid resolution for
#' interpolating repulsive forces (larger is slower but more accurate). A value
#' of zero disables interpolation, while a value of NULL only uses interpolation
#' for large datasets.
#' @param max_depth Integer scalar specifying the maximum depth of the
#' Barnes-Hut quad trees (larger is slower but more accurate). Set to a large
#' integer (e.g., 1000) to eliminate depth restrictions.
#' @param max_iter Integer scalar specifying the maximum number of iterations to
#' perform.
#' @inheritParams runUMAP
#' @inheritParams scran.chan::runTSNE.chan
#' @inherit runPCA return
#' @seealso [runTSNE.chan][scran.chan::runTSNE.chan]
#' @export
#' @rdname runTSNE
runTSNE.default <- function(object, perplexity = 30L, interpolate = NULL,
                            max_depth = 7L, max_iter = 500L,
                            ...,
                            approximate = TRUE, seed = 123456L,
                            threads = NULL) {
    assert_number(perplexity)
    assert_number(interpolate, null_ok = TRUE)
    assert_number(max_depth)
    assert_number(max_iter)
    assert_number(seed)
    scran.chan::runTSNE.chan(
        x = object,
        perplexity = perplexity,
        interpolate = interpolate,
        max.depth = max_depth,
        max.iter = max_iter,
        seed = as.integer(seed),
        approximate = approximate,
        downsample = NULL,
        num.threads = set_threads(threads)
    )
}
