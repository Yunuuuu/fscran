.get_mat_from_sce <- function(x, exprs_values, dimred, n_dimred) {
    if (!is.null(dimred)) {
        mat <- SingleCellExperiment::reducedDim(x, dimred)
        if (!is.null(n_dimred)) {
            if (length(n_dimred) == 1L) n_dimred <- seq_len(n_dimred)
            mat <- mat[, n_dimred, drop = FALSE]
        }
        mat
    } else {
        SummarizedExperiment::assay(x, exprs_values)
    }
}
