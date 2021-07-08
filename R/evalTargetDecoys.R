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
#' exampleFiles <- list.files(system.file('extdata', package = 'mzID'),
#'                           pattern = '*.mzid', full.names = TRUE)
#' mzIDexample <- mzID(exampleFiles[[2]])
#' decoyPlots<-evalTargetDecoys(object = mzIDexample, decoy = "isdecoy", score = "x!tandem:expect", log10 =T RUE, nBins = 30)
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
    # require some packages
    # require(mzID)
    # require(mzR)
    # require(dplyr)
    # require(ggplot2)
    # require(graphics)
    # require(ggpubr)
    # require(shiny)

    # check object class
    if (is(object, "mzID")) {
        df <- mzID::flatten(object)
    } else if (is.data.frame(object)) {
        df <- object
    } else if (is(object, "mzRident")) {
        df <- cbind(psms(object), score(object)[, -1])
    } else {
        "object should be of the class mzID, mzRident or dataframe"
    }

    # if one or more arguments in the function are missing, .select() is called, a pop-up window opens and the variables must be chosen manually.
    if (missing(decoy) | missing(score) | missing(log10)) {
        if (missing(log10)) log10 <- TRUE
        out <- .select(df, decoy, score, log10, nBins)
        decoy <- out$selDecoy # categorical
        score <- out$selScore # continu
        log10 <- out$log
        nBins <- out$nBins
    }

    # subset of dataframe
    data <- df[, c(decoy, score)]
    names(data) <- c("decoy", "score")
    data <- na.exclude(data)
    data$score <- as.double(data$score)

    # if variable 'score' is a character, change to continuous
    if (is.character(data$score)) {
        data$score <- as.numeric(as.character(data$score))
    }

    # check whether the selected variables are of the correct class
    if (!is.numeric(data$score)) stop("Score is not numeric")
    if (!is.logical(data$decoy)) stop("Decoy is not logical")

    # perform log10-transformation on variable 'score' if so indicated
    if (log10) {
        data$score <- -log10(as.numeric(data$score))
    }

    #############
    ### Plots ###
    #############

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


    out <- list(
        ppPlot = p1$ppPlot,
        histogram = p2,
        ppPlotZoom = p1$ppPlot + ylim(p1$yLim[1], p1$yLim[2]),
        histogramZoom = p2 + coord_cartesian(xlim = c(min(data$score), max(data$score[data$decoy])), expand = TRUE)
    )
    out$together <- ggarrange(plotlist = out, ncol = 2, nrow = 2)
    return(out)
}


## Create PP-plot
.ppPlot <- function(data, ylim = NULL) {
    pi0 <- sum(data$decoy) / sum(!data$decoy)

    x <- data$score[!data$decoy]
    Ft <- ecdf(x)
    Fd <- ecdf(data$score[data$decoy])
    df <- data_frame(Fdp = Fd(x), Ftp = Ft(x))
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
