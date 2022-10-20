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

## Helper to compute data for PP-plot
#' @importFrom stats ecdf
.ppData <- function(data, maxLength=1000) {
    pi0 <- sum(data$decoy) / sum(!data$decoy)

    x <- data$score[!data$decoy]
    Ft <- ecdf(x)
    Fd <- ecdf(data$score[data$decoy])
    if (length(x) > maxLength) {
      x <- quantile(x,seq(0,1,length=maxLength))
    }
    z <- Ft(x) - pi0 * Fd(x)
    df <- data.frame(Fdp = Fd(x), Ftp = Ft(x), z = z)
    #ylimHlp <- mean(Fd(x) == 1)
    ylimHlp <- Ft(max(data$score[data$decoy]))
    list(df = df, pi0 = pi0, ylimHlp = ylimHlp)
}


# Function for ggplot2-like colour scale in HCL space
gg_color_hue <- function(n) {
    hues <- seq(15, 375, length.out = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}
