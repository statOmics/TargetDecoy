#' Evaluate assumptions of the Target Decoys approach
#'
#' Create diagnostic plots to evaluate the TDA assumptions.
#' A histogram and PP plot allow to check both necessary assumptions.
#'
#' @param object A `data.frame`, \linkS4class{mzID} or \linkS4class{mzR} object.
#' @param decoy `logical`, to indicate if the peptide matches to a target or to
#'   a decoy.
#' @param score `numeric`, indicating the score of the peptide match, obtained
#'   by the search engine.
#' @param log10 `logical` to indicate if the score should be
#'   `-log10`-transformed.
#' @param nBins `numeric` indicating the number of bins in the histogram.
#'
#' @return
#' A list containing the PP-plot and histogram, a zoom of both plots and an
#' overview of all four plots.
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
#' decoyPlots <- evalTargetDecoys(mzIDexample,
#'     decoy = "isdecoy", score = "x\\!tandem:expect", log10 = TRUE
#' )
#'
#' # Plot only the overview of the four plots
#' decoyPlots$together
#'
#' # If you do not know the name of the score and/or decoy variable,
#' # or if you want to evaluate how many bins you want to use in the histogram,
#' # or if -log10 transformation is needed you can launch a shiny app by only
#' # specifying the mzID object
#'
#' # evalTargetDecoys(mzIDexample)
evalTargetDecoys <- function(object,
                             decoy = NULL,
                             score = NULL,
                             log10 = NULL,
                             nBins = 50) {
    ## Prepare score table
    data <- decoyScoreTable(
        object = object,
        decoy = decoy, score = score, log10 = log10
    )

    # create PP-plot
    p1 <- .ppPlot(data)

    # create histogram
    p2 <- ggplot(data, aes(score, fill = decoy, col = I("black"))) +
        geom_histogram(alpha = 0.5, bins = nBins, position = "identity") +
        labs(
            x = head(score), y = "",
            title = "Histogram of targets and decoys"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(size = rel(1.5)),
            axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            axis.title.y = element_text(angle = 0)
        )

    ppPlotZoom <- p1$ppPlot + ylim(p1$yLim[1], p1$yLim[2])
    histogramZoom <- p2 +
        coord_cartesian(
            xlim = c(min(data$score), max(data$score[data$decoy])),
            expand = TRUE
        )

    out <- list(
        ppPlot = p1$ppPlot,
        histogram = p2,
        ppPlotZoom = ppPlotZoom,
        histogramZoom = histogramZoom
    )
    out$together <- ggpubr::ggarrange(plotlist = out, ncol = 2, nrow = 2)
    return(out)
}


## Create PP-plot
.ppPlot <- function(data, ylim = NULL) {
    pi0 <- sum(data$decoy) / sum(!data$decoy)

    x <- data$score[!data$decoy]
    Ft <- ecdf(x)
    Fd <- ecdf(data$score[data$decoy])
    df <- data.frame(Fdp = Fd(x), Ftp = Ft(x))
    ylimHlp <- mean(Fd(x) == 1)

    ppPlot <- ggplot(df) +
        geom_point(aes(Fdp, Ftp), color = "dark gray") +
        geom_abline(slope = pi0, color = "black") +
        labs(title = "PP plot") +
        coord_cartesian(xlim = c(0, 1), ylim = ylim, expand = TRUE) +
        theme_bw() +
        theme(
            plot.title = element_text(size = rel(1.5)),
            axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            axis.title.y = element_text(angle = 0)
        )

   list(
        ppPlot = ppPlot,
        yLim = c(0, (1 - ylimHlp)),
        Fd = Fd,
        Ft = Ft
    )
}
