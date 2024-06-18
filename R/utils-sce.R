.get_mat_from_sce <- function(x, assay, dimred, n_dimred) {
    if (!is.null(dimred)) {
        # value is expected to be a matrix or matrix-like object with number of
        # rows equal to ncol(x).
        mat <- SingleCellExperiment::reducedDim(x, dimred)
        if (!is.null(n_dimred)) {
            if (length(n_dimred) == 1L) n_dimred <- seq_len(n_dimred)
            mat <- mat[, n_dimred, drop = FALSE]
        }
        # we always regard features in row and cells in column
        t(mat)
    } else {
        SummarizedExperiment::assay(x, assay)
    }
}
