#' Perform Principal components analysis with `scran.chan`
#'
#' @param object A matrix-like object. Rows are features and columns are cells.
#' @param ... Additional arguments passed to default methods.
#' @export
runPCA <- function(object, ...) UseMethod("runPCA")

#' @param assay Integer scalar or string indicating which assay of x
#' contains the expression values.
#' @param dimred String or integer scalar specifying the existing dimensionality
#' reduction results to use.
#' @param n_dimred Integer scalar or vector specifying the dimensions to use if
#' dimred is specified.
#' @param size_factors A numeric vector of length equal to the number of cells
#' in x, containing positive size factors for all cells.
#' @param name String specifying the name to be used to store the result in the
#' [reducedDims][SingleCellExperiment::reducedDim] or
#' [reductions][SeuratObject::Seurat-class] of the output.
#' @export
#' @rdname runPCA
runPCA.SingleCellExperiment <- function(object, ...,
                                        size_factors = NULL,
                                        assay = "logcounts",
                                        dimred = NULL, n_dimred = NULL,
                                        name = "PCA") {
    mat <- .get_mat_from_sce(object, assay, dimred, n_dimred)
    pca <- runPCA(object = mat, ...)
    SingleCellExperiment::reducedDim(object, name) <- pca
    object
}

#' @param layer Name of the layer to get from the assay data.
#' @export
#' @rdname runPCA
runPCA.Seurat <- function(object, ...,
                          assay = NULL, layer = "data",
                          dimred = NULL, n_dimred = NULL,
                          name = "PCA") {
    mat <- .get_mat_from_seurat(object, assay, layer, dimred, n_dimred)
    pca <- runPCA(object = mat, ...)
    reduction_key <- SeuratObject::Key(name, quiet = TRUE)
    rownames(pca) <- rownames(mat)
    colnames(pca) <- paste0(reduction_key, seq_len(ncol(pca)))
    object[[name]] <- SeuratObject::CreateDimReducObject(
        embeddings = pca,
        stdev = as.numeric(apply(pca, 2L, stats::sd)),
        assay = .get_assay_from_seurat(object, assay, layer, dimred),
        key = reduction_key
    )
    object
}

#' @param d Integer scalar specifying the number of top PCs to obtain.
#' @param scale Logical scalar indicating whether to scale rows to unit
#' variance.
#' @param subset_row Integer, logical or character vector specifying which
#' features to use in the PCA (e.g., highly variable genes). If NULL, all
#' features in x are used.
#' @param batch Vector or factor of length equal to the number of cells,
#' specifying the batch of origin for each cell. Alternatively NULL if all cells
#' belong to the same batch.
#' @param batch_method String indicating how batch should be handled (if it is
#' supplied). `"block"` is equivalent to linear regression on x prior to PCA,
#' while `"weight"` will only weight each batch so that they contribute equally
#' to the PCA.
#' @param force_integer Logical scalar indicating whether double-precision x
#' should be forced into integers.
#' @param no_sparse_copy Logical scalar indicating whether we should avoid a
#' copy when object is a [dgCMatrix][Matrix::dgCMatrix-class] This is more
#' memory efficient if the data has already been loaded into memory. If TRUE,
#' any setting of force.integer is ignored.
#' @param threads Integer scalar specifying the number of threads to use. If
#' `NULL`, all detected threads will be used. See
#' [detectCores][parallel::detectCores].
#' @seealso [runPCA.chan][scran.chan::runPCA.chan]
#' @return
#'  - `default` method: A numeric matrix where rows are cells and columns are
#'                    the two dimensions of the embedding.
#'  - `SingleCellExperiment` method: embedding was added into
#'    [reducedDims][SingleCellExperiment::reducedDim] named as `name`.
#'  - `Seurat` method: embedding was added into
#'    [reductions][SeuratObject::Seurat-class] named as `name`.
#' @export
#' @rdname runPCA
runPCA.default <- function(object, d = 50L, scale = FALSE, ...,
                           subset_row = NULL, batch = NULL, batch_method = NULL,
                           force_integer = TRUE, no_sparse_copy = TRUE,
                           threads = NULL) {
    threads <- set_threads(threads)

    # dimensionality reduction -----------------------------------------
    batch_pcs <- scran.chan::runPCA.chan(
        x = scran.chan::initializeSparseMatrix(
            object,
            force.integer = force_integer,
            no.sparse.copy = no_sparse_copy,
            by.column = TRUE,
            num.threads = threads
        ),
        num.comp = as.integer(d),
        subset = subset_row,
        scale = scale,
        num.threads = threads,
        batch = batch,
        batch.method = batch_method,
        rotation = TRUE
    )
    # batch_pcs$components:
    # a numeric matrix containing the top principal components. Each row
    # corresponds to a PC and each column corresponds to a cell.
    out <- t(batch_pcs$components)
    attr(out, "percentVar") <- batch_pcs$prop.variance
    attr(out, "rotation") <- batch_pcs$rotation
    out
}
