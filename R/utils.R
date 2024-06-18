set_threads <- function(threads,
                        arg = rlang::caller_arg(threads),
                        call = rlang::caller_call()) {
    assert_(threads,
        function(x) is_scalar_numeric(x) && x > 0L,
        "a positive number",
        null_ok = TRUE, arg = arg, call = call
    )
    if (is.null(threads)) {
        parallel::detectCores()
    } else {
        as.integer(threads)
    }
}

is_scalar_numeric <- function(x) length(x) == 1L && is.numeric(x)
