#' Score marker genes
#' @inheritParams runPCA
#' @export
scoreMarkers <- function(object, ...) UseMethod("scoreMarkers")

#' @export
#' @rdname scoreMarkers
scoreMarkers.SingleCellExperiment <- function(object, ..., assay = "logcounts") {
    mat <- .get_mat_from_sce(object, assay)
    scoreMarkers(object = mat, ...)
}

#' @export
#' @rdname scoreMarkers
scoreMarkers.Seurat <- function(object, ..., assay = NULL, layer = "data") {
    mat <- .get_mat_from_seurat(object, assay, layer)
    scoreMarkers(object = mat, ...)
}

#' @param groups A vector specifying the group assignment for each cell in x.
#' @inheritParams scran.chan::scoreMarkers.chan
#' @param simple_means_only Logical scalar indicating whether to only report the
#' means for the simple effect sizes, i.e., log-fold change and delta-detected.
#' @param sort_by String specifying the column to use for sorting genes in
#' descending order (except if it ends with `.rank`, in which case it is sorted
#' in ascending order). If `NULL`, no sorting is performed.
#' @param all_pairwise Logical scalar indicating whether to report the full
#' effects for every pairwise comparison.
#' @inheritParams logNormCounts
#' @inherit scran.chan::scoreMarkers.chan details return
#' @rdname scoreMarkers
#' @export
scoreMarkers.default <- function(object, groups, lfc = 0L,
                                 simple_means_only = TRUE,
                                 sort_by = "cohen.rank", all_pairwise = FALSE,
                                 ...,
                                 batch = NULL,
                                 force_integer = TRUE,
                                 no_sparse_copy = TRUE,
                                 threads = NULL) {
    rlang::check_dots_empty()
    scran.chan::scoreMarkers.chan(
        x = scran.chan::initializeSparseMatrix(
            object,
            force.integer = force_integer,
            no.sparse.copy = no_sparse_copy,
            by.column = TRUE,
            num.threads = threads
        ),
        groups = groups, batch = batch,
        lfc = lfc,
        num.threads = set_threads(threads),
        simple.means.only = simple_means_only,
        sort.by = sort_by,
        all.pairwise = all_pairwise
    )
}
