#########################################################################################
#########################################################################################
# Function to create all the PP plots in one figure for scores from files from one engine
#########################################################################################
#########################################################################################
# mzObjects: list of mzID files
# decoy: string,  name of the variable isdecoy in this mzID file (typical "isdecoy" or "isDecoy")
# score: string, name of the variabele score in this mzID file
# log10: boolean, TRUE if the scores from the mzID are returned log10 transformed
# sameSearchEngine: boolean, TRUE if scores result from the same search engine

#' Create all the PP plots in one figure for scores from multiple objects
#'
#' @param object_list
#' @inheritParams decoyScoreTable
#' @param sameSearchEngine
#'
#' @return
#' @export
#'
#' @examples
#' library(mzID)
#'
#' ## Use two example files from the mzID package
#' exampleFiles <- system.file(
#'     "extdata", c("55merge_omssa.mzid", "55merge_tandem.mzid"),
#'     package = "mzID"
#' )
#' mzObjects <- lapply(exampleFiles, mzID)
#'
#' createPPlotObjects(mzObjects,
#'     decoy = "isdecoy",
#'     score = c("omssa:evalue", "x\\!tandem:expect"),
#'     log10 = TRUE
#' )
createPPlotObjects <- function(object_list, decoy = NULL, score = NULL, log10 = TRUE, sameSearchEngine = TRUE) {
    if (!is.list(object_list)) object_list <- list(object_list)

    # TODO: do we need this? The Shiny app is already called through processObjects?
    #       Also, is `sameSearchEngine` really necessary?
    # if ((is.null(decoy) | is.null(score) | is.null(log10)) & sameSearchEngine) {
    #     object <- mzObjects[[1]]
    #     df <- .getDF(object)
    #     if (missing(log10)) log10 <- TRUE
    #     out <- .select(df, decoy, score, log10)
    #     decoy <- out$selDecoy # categorical
    #     score <- out$selScore # continu
    #     log10 <- out$log
    # }

    decoyScoreTables <- processObjects(object_list, decoy, score, log10)

    if (is.null(names(decoyScoreTables))) {
        names(decoyScoreTables) <- 1:length(decoyScoreTables)
    }

    ppData <- lapply(
        1:length(decoyScoreTables),
        function(i, tables) {
            pi0 <- sum(tables[[i]]$decoy) / sum(!tables[[i]]$decoy)
            x <- tables[[i]]$score[!tables[[i]]$decoy]
            Ft <- ecdf(x)
            Fd <- ecdf(tables[[i]]$score[tables[[i]]$decoy])
            df <- tibble(Fdp = Fd(x), Ftp = Ft(x), z = Ft(x) - pi0 * Fd(x), pi0 = pi0, search = names(tables)[i])
            return(df)
        },
        tables = decoyScoreTables
    ) %>% bind_rows()

    # Plot with original pi0 line and points
    plot1 <- ppData %>%
        ggplot(aes(Fdp, Ftp, color = search)) +
        geom_point() +
        ggtitle("PP plot") +
        ylab("FTarget") +
        xlab("FDecoy") +
        theme_bw() +
        geom_abline(
            slope = ppData %>% pull(pi0) %>% unique(),
            color = gg_color_hue(length(object_list))
        )

    # Plot with rescaled pi0 line on 0 and points
    plot2 <- ppData %>%
        ggplot(aes(Fdp, z, color = search)) +
        geom_point() +
        ggtitle("PP plot - pi0 subtracted") +
        ylab("FTarget-pi0") +
        xlab("FDecoy") +
        theme_bw() +
        geom_abline(slope = 0)
    return(list(plot1, plot2))
}


# Helper to make decoy score tables from multiple mzID objects
processObjects <- function(object_list, decoy, score, log10) {
    arg_list <- .check_args(
        object_list = object_list,
        decoy = decoy, score = score, log10 = log10
    )
    object_list <- arg_list$object_list
    decoy <- arg_list$decoy
    score <- arg_list$score
    log10 <- arg_list$log10

    out <- vector("list", length = length(object_list))
    for (i in seq_along(object_list)) {
        out[[i]] <- decoyScoreTable(
            object = object_list[[i]],
            decoy = decoy[[i]],
            score = score[[i]],
            log10 = log10[[i]]
        )
    }
    names(out) <- names(object_list)
    out
}
