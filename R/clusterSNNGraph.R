#' Graph-based clustering with `scran.chan`
#'
#' @seealso [clusterSNNGraph.chan][scran.chan::clusterSNNGraph.chan]
#' @export
clusterSNNGraph <- function(object, ...) UseMethod("clusterSNNGraph")

#' @inheritParams runPCA
#' @export
#' @rdname clusterSNNGraph
clusterSNNGraph.SingleCellExperiment <- function(object, ...,
                                                 dimred = "PCA",
                                                 n_dimred = NULL,
                                                 assay = NULL) {
    mat <- .get_mat_from_sce(object, assay, dimred, n_dimred)
    clusterSNNGraph(object = mat, ...)
}

#' @inheritParams runPCA
#' @export
#' @rdname clusterSNNGraph
clusterSNNGraph.Seurat <- function(object, ...,
                                   dimred = "PCA", n_dimred = NULL,
                                   assay = NULL, layer = NULL) {
    mat <- .get_mat_from_seurat(object, assay, layer, dimred, n_dimred)
    clusterSNNGraph(object = mat, ...)
}

#' @inheritParams runMNN
#' @param method String specifying the community detection method to use.
#' Options are multi-level (`"multilevel"`), Walktrap (`"walktrap"`) or Leiden
#' (`"leiden"`).
#' @param scheme String specifying the weighting scheme to use for constructing
#' the SNN graph. This can be `"rank"` (default), `"jaccard"` or `"number"`.
#' @param resolution Numeric scalar specifying the resolution to use for
#' multi-level or Leiden clustering.
#' @param steps Integer scalar specifying the number of steps to use for
#' Walktrap clustering.
#' @param objective String specifying the objective function to use for Leiden
#' clustering: "CPM" or "modularity".
#' @inheritParams scran.chan::clusterSNNGraph.chan
#' @param seed Integer scalar specifying the seed to use for multi-level or
#' Leiden clustering.
#' @return An integer vector with cluster assignments for each cell. Each method
#' may also return additional attributes.
#' - For method=`"multilevel"`, we have:
#'    * `levels`, a list of integer vectors with cluster assignments for each
#' cell at each level. Assignments are sorted by decreasing resolution (i.e.,
#' fewer, larger clusters).
#'    * `modularity`, a numeric vector containing the modularity of each level.
#'    * `best`, the level with the lowest modularity.
#' - For method=`"leiden"`, we have:
#'    * `quality`: a numeric scalar containing the quality of the clustering
#'      (either the modularity or a related score).
#' - For method=`"walktrap"`, we have:
#'    * `merges`: an integer matrix specifying how the clusters were merged to
#'      obtain membership. Each row corresponds to a merge step and contains the
#'      IDs of the temporary clusters (not the same as those in membership).
#'    * `modularity`: a numeric vector containing the modularity before and
#'      after each merge step.
#' @export
#' @rdname clusterSNNGraph
clusterSNNGraph.default <- function(object, k = 15L, method = "leiden", ...,
                                    scheme = NULL,
                                    resolution = 1L,
                                    objective = "modularity", steps = 4L,
                                    approximate = TRUE, seed = 123456L,
                                    threads = 1L) {
    assert_number(k)
    method <- match.arg(method, c("multilevel", "walktrap", "leiden"))
    assert_number(resolution)
    assert_number(steps)
    objective <- match.arg(objective, c("CPM", "modularity"))
    assert_bool(approximate)
    assert_number(seed)
    # run MNN --------------------------------------------------------
    clustering <- scran.chan::clusterSNNGraph.chan(
        x = object,
        num.neighbors = as.integer(k), ...,
        weight.scheme = scheme, method = method,
        resolution = resolution, objective = objective,
        approximate = approximate,
        num.threads = as.integer(threads),
        downsample = NULL,
        drop = TRUE, seed = as.integer(seed)
    )
    out <- clustering$membership
    for (i in setdiff(names(clustering), "membership")) {
        attr(out, i) <- clustering[[i]]
    }
    out
}
