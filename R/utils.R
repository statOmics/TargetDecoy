## Check that provided arguments have same length, expand any that have length 1
.check_args <- function(...) {
    arg_list <- list(...)
    arg_lengths <- lengths(arg_list)

    multi_args <- which(arg_lengths > 1)
    if (any(multi_args)) {
        multi_arg_lengths <- arg_lengths[multi_args]

        if (length(unique(multi_arg_lengths)) > 1) {
            bad_args <- paste(names(arg_list[multi_args]), collapse = ", ")
            stop("Arguments `", bad_args, "` have different lenghts.",
                "\n  They should have the same length or length == 1.",
                call. = FALSE
            )
        }

        ## Args with length > 1 guaranteed to have same length at this point
        required_length <- unique(multi_arg_lengths)

        ## Expand arguments with length 1
        single_args <- which(arg_lengths == 1)
        arg_list[single_args] <- lapply(arg_list[single_args],
            rep.int,
            times = required_length
        )
    }

    arg_list
}
