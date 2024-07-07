initialize_matrix <- function(x, ...) UseMethod("initialize_matrix")

#' @export
initialize_matrix.default <- function(x, ..., by_column = TRUE,
                                      force_integer = TRUE,
                                      no_sparse_copy = TRUE, threads = NULL) {
    scran.chan::initializeSparseMatrix(
        x,
        force.integer = force_integer,
        no.sparse.copy = no_sparse_copy,
        by.column = by_column,
        num.threads = set_threads(threads)
    )
}

#' @export
initialize_matrix.list <- function(x, ...) {
    if (!setequal(names(x), c("pointer", "rownames", "colnames"))) {
        cli::cli_abort(c(
            "invalid {.arg x}",
            i = paste(
                "{.arg x} must be a list with 3 elements:",
                "pointer, rownames, and colnames"
            )
        ))
    } else if (!inherits(.subset2(x, "pointer"), "externalptr")) {
        cli::cli_abort("{.code $pointer} must be a {.cls externalptr}")
    } else if (!is.null(.subset2(x, "rownames")) &&
        !is.character(.subset2(x, "rownames"))) {
        cli::cli_abort("{.code $rownames} must be a {.cls character} or `NULL`")
    } else if (!is.null(.subset2(x, "colnames")) &&
        !is.character(.subset2(x, "colnames"))) {
        cli::cli_abort("{.code $colnames} must be a {.cls character} or `NULL`")
    }
    x
}
