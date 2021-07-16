#' Evaluate assumptions of the Target Decoys Approach for multiple search
#' engines
#'
#' Create diagnostic PP plots in one figure to evaluate the TDA assumptions for
#' multiple search engines. The function provides the possibility to evaluate
#' each of the sub-engines and the overall itself.
#'
#' @inheritParams decoyScoreTable
#' @param scores A `character` vector of multiple score names from the input
#'     object. Typically from different search engines.
#'
#' @return
#' One PP plot with all original pi0, and a standardized / rescaled PP plot with
#' all `pi0` set to 0.
#'
#' @author Elke Debrie, Lieven Clement
#'
#' @export
#'
#' @examples
#' library(mzID)
#'
#' ## Use one of the example files in the mzID package
#' exampleFile <- system.file("extdata", "55merge_tandem.mzid", package = "mzID")
#' mzIDexample <- mzID(exampleFile)
#'
#' plots <- createPPlotScores(mzIDexample,
#'     scores = c("x\\!tandem:hyperscore", "x\\!tandem:expect"),
#'     decoy = "isdecoy", log10 = TRUE
#' )
createPPlotScores <- function(object, scores, decoy = NULL, log10 = TRUE) {
    tables <- processScores(
        object = object,
        scores = scores, decoy = decoy, log10 = log10
    )
    ppData <- ppScoresData(tables)
    ppScoresPlots(ppData)
}


# Helper to make decoy score tables from a single mzID object for multiple
# scores.
processScores <- function(object, scores, decoy, log10) {
    arg_list <- .check_args(scores = scores, log10 = log10)
    scores <- arg_list$scores
    log10 <- arg_list$log10

    out <- vector("list", length = length(scores))
    for (i in seq_along(scores)) {
        out[[i]] <- decoyScoreTable(
            object = object,
            decoy = decoy,
            score = scores[[i]],
            log10 = log10[[i]]
        )
    }
    names(out) <- scores
    out
}


ppScoresData <- function(tables) {
    tmp <- lapply(tables, .ppData)
    dfs <- lapply(tmp, `[[`, "df")
    df <- do.call(rbind, dfs)
    # Add `id` column from table names
    df$id <- rep(names(tables), vapply(dfs, nrow, integer(1)))
    if (is.null(df$id)) {
        df$id <- rep(seq_along(tables), vapply(dfs, nrow, integer(1)))
    }
    pi0 <- vapply(tmp, `[[`, FUN.VALUE = double(1), "pi0")
    list(df = df, pi0 = pi0)
}


ppScoresPlots <- function(ppData) {
    df <- ppData$df
    df$id <- as.factor(df$id)
    pi0 <- ppData$pi0
    base_plot <- ggplot(df) +
        xlab("FDecoy") +
        theme_bw()

    p1 <- base_plot +
        geom_point(aes(Fdp, Ftp, color = id)) +
        ggtitle("PP plot") +
        ylab("FTarget") +
        geom_abline(
            slope = pi0,
            color = gg_color_hue(length(pi0))
        )

    p2 <- base_plot +
        geom_point(aes(Fdp, z, color = id)) +
        ggtitle("PP plot - pi0 subtracted") +
        ylab("FTarget-pi0") +
        geom_abline(slope = 0)

    list(p1, p2)
}
