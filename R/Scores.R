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
#' exampleFile <- system.file('extdata', '55merge_tandem.mzid', package = 'mzID')
#' mzIDexample <- mzID(exampleFile)
#'
#' plots <- createPPlotScores(mzIDexample,
#'     scores = c("x\\!tandem:hyperscore", "x\\!tandem:expect"),
#'     decoy = "isdecoy", log10 = TRUE
#' )
createPPlotScores <- function(object, scores, decoy = NULL, log10 = TRUE) {
    decoyScoreTables <- processScores(
        object = object,
        scores = scores, decoy = decoy, log10 = log10
    )

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
    ) %>%
        bind_rows() %>%
        arrange(search)

    plot1 <- ppData %>%
        ggplot(aes(Fdp, Ftp, color = search)) +
        geom_point() +
        ggtitle("PP plot") +
        ylab("FTarget") +
        xlab("FDecoy") +
        theme_bw() +
        geom_abline(
            slope = ppData %>% pull(pi0) %>% unique(),
            color = gg_color_hue(length(scores))
        )

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


# Function to make a list of the processed mzIDFiles
# returns a list of tables with decoys (TRUE/FALSE) and scores from a given mzID file

processScores <- function(object, scores, decoy = NULL, log10 = TRUE) {
    if (!is.null(decoy) & length(decoy) == 1) decoy <- rep(decoy, length(scores))
    if (!is.null(log10) & length(log10) == 1) log10 <- rep(log10, length(scores))

    decoyScoreTables <- lapply(
        1:length(scores),
        function(i, decoy, scores, log10, object) {
              decoyScoreTable(object, decoy[i], scores[i], log10[i])
          },
        decoy = decoy,
        scores = scores,
        log10 = log10,
        object = object
    )
    names(decoyScoreTables) <- scores
    return(decoyScoreTables)
}

# Function for ggplot2-like colour scale in HCL space

gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
