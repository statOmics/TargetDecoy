#' Prepare score table for decoys
#'
#' Takes an input object and returns a score table for the decoys.
#'
#' @param object A `data.frame`, \linkS4class{mzID} or \linkS4class{mzRident} object.
#' @param decoy `character`, name of the variable that indicates if the peptide
#'   matches to a target or to a decoy.
#' @param score `numeric`, indicating the score of the peptide match, obtained
#'   by the search engine.
#' @param log10 `logical` to indicate if the score should be
#'   `-log10`-transformed.
#' @param nBins `numeric` indicating the number of bins in the histogram.
#'
#' @return
#' A `data.frame` with a logical `"decoy"` column and numeric `"scores"`.
#'
#' @author Elke Debrie, Lieven Clement
#' @keywords internal
decoyScoreTable <- function(object, decoy, score, log10 = TRUE) {
    df <- .getDF(object)

    is_missing <- !(c(decoy, score) %in% colnames(df))
    if (any(is_missing)) {
        missing <- paste(c(decoy, score)[is_missing], collapse = ", ")
        stop("Variable(s) not found: ", missing, call. = FALSE)
    }

    table <- df[, c(decoy, score)]
    names(table) <- c("decoy", "score")
    table <- stats::na.exclude(table)
    table$score <- as.double(table$score)

    # if variable 'score' is a character, change to continuous
    if (is.character(table$score)) {
        table$score <- as.numeric(as.character(table$score))
    }

    if (!is.logical(table$decoy)) stop("`decoy` is not logical.", call. = FALSE)

    # perform log10-transformation on variable 'score' if so indicated
    if (log10) {
        table$score <- -log10(as.numeric(table$score))
    }

    return(table)
}


#' @importFrom mzID flatten
#' @importFrom mzR psms score
.getDF <- function(object) {
    # check object class
    if (is.data.frame(object)) {
        return(object)
    } else if (is(object, "mzID")) {
        df <- flatten(object)
    } else if (is(object, "mzRident")) {
        df <- cbind(psms(object), score(object)[, -1])
    } else {
        stop(
            "`object` should be of class 'mzID', 'mzRident' or 'data.frame',",
            "\n  not '", class(object), "'.", call. = FALSE
        )
    }
    df
}


#' Evaluate assumptions of the Target Decoys approach
#'
#' Create diagnostic plots to evaluate the TDA assumptions.
#' A histogram and PP plot allow to check both necessary assumptions.
#'
#' @inheritParams decoyScoreTable
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
#'
#'
#' ## mzRident objects can also be used
#' library(mzR)
#'
#' if (requireNamespace("msdata", quietly = TRUE)) {
#'    ## Using example file from msdata
#'    file <- system.file("mzid", "Tandem.mzid.gz", package="msdata")
#'    mzid <- openIDfile(file)
#' }
#' decoyPlots2 <- evalTargetDecoys(mzid,
#'     decoy = "isDecoy", score = "X.Tandem.expect", log10 = TRUE
#' )
#' decoyPlots2$together
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
            x = utils::head(score), y = "",
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



evalTargetDecoysPPPlot <- function(object,
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

    out <- ppPlot = p1$ppPlot

    return(out)
}

evalTargetDecoysHist <- function(object,
                             decoy = NULL,
                             score = NULL,
                             log10 = NULL,
                             nBins = 50) {
    ## Prepare score table
    data <- decoyScoreTable(
        object = object,
        decoy = decoy, score = score, log10 = log10
    )

    # create histogram
    p2 <- ggplot(data, aes(score, fill = decoy, col = I("black"))) +
        geom_histogram(alpha = 0.5, bins = nBins, position = "identity") +
        labs(
            x = utils::head(score), y = "",
            title = "Histogram of targets and decoys"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(size = rel(1.5)),
            axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            axis.title.y = element_text(angle = 0)
        )



    out <- histogram = p2
    return(out)
}




evalTargetDecoysZoomPP <- function(object,
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


    ppPlotZoom <- p1$ppPlot + ylim(p1$yLim[1], p1$yLim[2])

    out <- ppPlotZoom = ppPlotZoom
    return(out)
}

evalTargetDecoysZooomHist <- function(object,
                             decoy = NULL,
                             score = NULL,
                             log10 = NULL,
                             nBins = 50) {
    ## Prepare score table
    data <- decoyScoreTable(
        object = object,
        decoy = decoy, score = score, log10 = log10
    )


    # create histogram
    p2 <- ggplot(data, aes(score, fill = decoy, col = I("black"))) +
        geom_histogram(alpha = 0.5, bins = nBins, position = "identity") +
        labs(
            x = utils::head(score), y = "",
            title = "Histogram of targets and decoys"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(size = rel(1.5)),
            axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            axis.title.y = element_text(angle = 0)
        )

    histogramZoom <- p2 +
        coord_cartesian(
            xlim = c(min(data$score), max(data$score[data$decoy])),
            expand = TRUE
        )

    out <- histogramZoom = histogramZoom

    return(out)
}


evalTargetDecoysTogether <- function(object,
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
            x = utils::head(score), y = "",
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
    return(out$together)
}
