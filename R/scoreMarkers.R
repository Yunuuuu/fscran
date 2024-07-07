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
#' @return A list of data frame of marker statistics.
#' Each data frame corresponds to a group in `groups` and contains:
#'
#' - `mean`: the mean expression across all cells in the current group.
#' - `detected`: proportion of cells with detectable expression in the current
#'             group.
#' - `logFC`: the mean of the log-fold changes in expression compared to other
#'          groups.
#' - `delta.detected`: the mean of the difference in the detected proportions
#'                   compared to other groups.
#' - `cohen.min`: the smallest Cohen's d across all pairwise comparisons
#'              involving the current group.
#' - `cohen.mean`: the mean Cohen's d across all pairwise comparisons involving
#'               the current group.
#' - `cohen.rank`: the minimum rank of the Cohen's d across all pairwise
#'               comparisons.
#' - `auc.min`: the smallest AUC across all pairwise comparisons involving the
#'              current group.
#' - `auc.mean`: the mean AUC across all pairwise comparisons involving the
#'               current group.
#' - `auc.rank`: the minimum rank of the AUC across all pairwise comparisons.
#'
#' Rows are sorted by the specified column in `sort.by`.
#'
#' If `simple_means_only=FALSE`, `logFC` and `delta.detected` are
#' renamed to `logFC.mean` and `delta.detected.mean`, respectively.
#' In addition, the corresponding `*.min` and `*.rank` columns are
#' also reported.
#' These can be interpreted in a similar manner as `cohen.min`,
#' `cohen.rank`, etc. for these effect sizes.
#'
#' If `batch` is supplied, this list will also contain `per.batch`.
#' This is a list containing `mean` and `detected`, each of which are
#' lists of data frames containing the batch-specific statistics for each group.
#' (In this case, `statistics` contains the averaged statistics for each
#' gene across batches)
#'
#' If `all_pairwise=TRUE`, this list will also contain `pairwise`, a list of
#' lists of data frames. Each data frame contains the statistics for the
#' pairwise comparison between groups, e.g., `pairwise$A$B` contains the
#' statistics for `A versus B` where large effects correspond to upregulation in
#' A. Note that rows correspond to the order in `object`. `sort_by` has no
#' effect on these data frames.
#'
#' @seealso [scoreMarkers.chan][scran.chan::scoreMarkers.chan]
#' @rdname scoreMarkers
#' @export
scoreMarkers.default <- function(object, groups, lfc = 0L,
                                 simple_means_only = TRUE,
                                 sort_by = "cohen.rank", all_pairwise = FALSE,
                                 ...,
                                 batch = NULL,
                                 force_integer = FALSE,
                                 no_sparse_copy = TRUE,
                                 threads = NULL) {
    rlang::check_dots_empty()
    threads <- set_threads(threads)
    out <- scran.chan::scoreMarkers.chan(
        x = initialize_matrix(
            object,
            force_integer = force_integer,
            no_sparse_copy = no_sparse_copy,
            threads = threads
        ),
        groups = groups, batch = batch,
        lfc = lfc,
        num.threads = threads,
        simple.means.only = simple_means_only,
        sort.by = sort_by,
        all.pairwise = all_pairwise
    )
    if (all_pairwise) out else out$statistics
}
