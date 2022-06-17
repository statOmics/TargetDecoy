#' Evaluate assumptions of the Target Decoys approach
#'
#' Create diagnostic plots to evaluate the TDA assumptions.
#' A histogram and PP plot allow to check both necessary assumptions.
#'
#' @inheritParams decoyScoreTable
#' @param nBins `numeric` indicating the number of bins in the histogram.
#'     When this value is missing, a Shiny gadget is launched to select it
#'     interactively.
#' @param zoom Logical value indicating whether a zoomed version of the plot
#'     should be returned. Default: `FALSE`.
#'
#' @return
#' `evalTargetDecoys` returns an overview of the following four plots:
#'
#'   1. A PP plot showing the empirical cumulative distribution of the target
#'     distribution in function of that of the decoy distribution
#'   2. A histogram showing the score distributions of the decoys and non-decoys
#'   3. A zoomed PP plot
#'   4. A zoomed histogram
#'
#' `evalTargetDecoysPPPlot` generates the PP plot only (1.) or the zoomed
#'  version (3.) if `zoom = TRUE`.
#'
#' `evalTargetDecoysHist` generates the histogram only (2.) or the zoomed
#'  version (4.) if `zoom = TRUE`.

#' @section The Shiny gadget:
#'
#' Sometimes the variable names are not known up front. If this is the case, the
#' `evalTargetDecoys*()` functions can be called with only an input object. This
#' launches a Shiny gadget that allows selecting the variables interactively. A
#' histogram and PP-plot of the selected variables are created on the fly for
#' previewing, together with a snapshot of the selected data.
#'
#' @author Elke Debrie, Lieven Clement, Milan Malfait
#'
#' @examples
#' library(mzID)
#'
#' ## Use one of the example files in the mzID package
#' exampleFile <- system.file("extdata", "55merge_tandem.mzid", package = "mzID")
#' mzIDexample <- mzID(exampleFile)
#'
#' # Plot the overview of the four plots
#' evalTargetDecoys(mzIDexample,
#'     decoy = "isdecoy", score = "x\\!tandem:expect", log10 = TRUE
#' )
#'
#' # Plot the PP plot only
#' evalTargetDecoysPPPlot(mzIDexample,
#'     decoy = "isdecoy", score = "x\\!tandem:expect", log10 = TRUE
#' )
#'
#' # Plot the zoomed PP plot only
#' evalTargetDecoysPPPlot(mzIDexample,
#'     decoy = "isdecoy", score = "x\\!tandem:expect", log10 = TRUE,
#'     zoom = TRUE
#' )
#'
#' # Plot the histogram only
#' evalTargetDecoysHist(mzIDexample,
#'     decoy = "isdecoy", score = "x\\!tandem:expect", log10 = TRUE
#' )
#'
#' # Plot the zoomed histogram only
#' evalTargetDecoysHist(mzIDexample,
#'     decoy = "isdecoy", score = "x\\!tandem:expect", log10 = TRUE,
#'     zoom = TRUE
#' )
#'
#' ## mzRident objects can also be used
#' library(mzR)
#'
#' if (requireNamespace("msdata", quietly = TRUE)) {
#'     ## Using example file from msdata
#'     file <- system.file("mzid", "Tandem.mzid.gz", package = "msdata")
#'     mzid <- openIDfile(file)
#' }
#' evalTargetDecoys(mzid,
#'     decoy = "isDecoy", score = "X.Tandem.expect", log10 = TRUE
#' )
#' @name evalTargetDecoys
NULL


#' @export
#' @rdname evalTargetDecoys
evalTargetDecoys <- function(object,
                             decoy = NULL, score = NULL,
                             log10 = TRUE, nBins = 50) {

    vars <- list(decoy = decoy, score = score, log = log10, nBins = nBins)
    if (any(vapply(vars, is.null, logical(1)))) {
        vars <- .select_vars(object,
            decoy = decoy, score = score,
            log10 = log10, nBins = nBins
        )
        decoy <- vars$decoy
        score <- vars$score
        log10 <- vars$log
        nBins <- vars$nBins
    }

    # create PP-plot
    ppPlot <- evalTargetDecoysPPPlot(
        object = object,
        decoy = decoy, score = score,
        log10 = log10, zoom = FALSE
    )

    # create histogram
    ppHist <- evalTargetDecoysHist(
        object = object,
        decoy = decoy, score = score,
        log10 = log10, nBins = nBins,
        zoom = FALSE
    )

    # create zoomed PP-plot
    ppPlotZoom <- evalTargetDecoysPPPlot(
        object = object,
        decoy = decoy, score = score,
        log10 = log10, zoom = TRUE
    )

    # create zoomed histogram
    histogramZoom <- evalTargetDecoysHist(
        object = object,
        decoy = decoy, score = score,
        log10 = log10, nBins = nBins,
        zoom = TRUE
    )

    ggpubr::ggarrange(
        ppPlot, ppHist, ppPlotZoom, histogramZoom,
        ncol = 2, nrow = 2
    )
}

#' @rdname evalTargetDecoys
#' @export
evalTargetDecoysPPPlot <- function(object,
                                   decoy = NULL, score = NULL,
                                   log10 = TRUE, zoom = FALSE) {

    vars <- list(decoy = decoy, score = score, log = log10)
    if (any(vapply(vars, is.null, logical(1)))) {
        vars <- .select_vars(object,
            decoy = decoy, score = score, log10 = log10
        )
        decoy <- vars$decoy
        score <- vars$score
        log10 <- vars$log
    }

    ## Prepare score table
    data <- decoyScoreTable(
        object = object,
        decoy = decoy, score = score, log10 = log10
    )

    # create PP-plot
    p <- .ppPlot(data)

    if (zoom) {
        out <- p$ppPlot + ylim(p$yLim[1], p$yLim[2])
    } else {
        out <- p$ppPlot
    }
    return(out)
}

## Helper to create PP-plot
.ppPlot <- function(data, ylim = NULL) {
    ppData <- .ppData(data)
    df <- ppData$df
    pi0 <- ppData$pi0
    ylimHlp <- ppData$ylimHlp

    ## Avoid R CMD check "no visible binding" warnings
    Fdp <- Ftp <- NULL

    ppPlot <- ggplot(df) +
        geom_point(aes(Fdp, Ftp), color = "dark gray", na.rm = TRUE) +
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

   list(ppPlot = ppPlot, yLim = c(0, (1 - ylimHlp)))
}


#' @rdname evalTargetDecoys
#' @export
evalTargetDecoysHist <- function(object,
                                 decoy = NULL, score = NULL,
                                 log10 = TRUE, nBins = 50,
                                 zoom = FALSE) {

    vars <- list(decoy = decoy, score = score, log = log10, nBins = nBins)
    if (any(vapply(vars, is.null, logical(1)))) {
        vars <- .select_vars(object,
            decoy = decoy, score = score,
            log10 = log10, nBins = nBins
        )
        decoy <- vars$decoy
        score <- vars$score
        log10 <- vars$log
        nBins <- vars$nBins
    }

    ## Prepare score table
    data <- decoyScoreTable(
        object = object,
        decoy = decoy, score = score, log10 = log10
    )

    # create histogram
    out <- ggplot(data, aes(score, fill = decoy, col = I("black"))) +
        geom_histogram(alpha = 0.5, bins = nBins, position = "identity") +
        labs(
            x = score, y = "",
            title = "Histogram of targets and decoys"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(size = rel(1.5)),
            axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            axis.title.y = element_text(angle = 0)
        )

    if (zoom) {
        out <- out +
            coord_cartesian(
                xlim = c(min(data$score), max(data$score[data$decoy])),
                expand = TRUE
            )
    }

    return(out)
 }
