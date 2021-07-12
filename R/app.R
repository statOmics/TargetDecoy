#' Prepare score table for decoys
#'
#' Takes an input object and returns a score table for the decoys. If one of the
#' arguments besides `object` is missing, will open a Shiny app to interactively
#' select the variables.
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
#' A `data.frame` with a logical `"decoy"` column and numeric `"scores"`.
#'
#' @author Elke Debrie, Lieven Clement
#' @keywords internal
#'
#' @examples
#' library(mzID)
#'
#' ## Use one of the example files in the mzID package
#' exampleFile <- system.file('extdata', '55merge_tandem.mzid', package = 'mzID')
#' mzIDexample <- mzID(exampleFile)
#' score_table <- decoyScoreTable(mzIDexample,
#'     decoy = "isdecoy", score = "x\\!tandem:expect", log10 = TRUE
#' )
decoyScoreTable <- function(object, decoy = NULL, score = NULL, log10 = TRUE) {
    df <- .getDF(object)

    if (is.null(decoy) || is.null(score) || is.null(log10)) {
        if (missing(log10)) log10 <- TRUE
        out <- .select(df, decoy, score, log10)
        decoy <- out$selDecoy # categorical
        score <- out$selScore # continu
        log10 <- out$log
    }

    if (!(decoy %in% colnames(df))) {
       stop("`decoy = '", decoy, "'` not found in input object.", call. = FALSE)
    }
    if (!(score %in% colnames(df))) {
       stop("`score = '", score, "'` not found in input object.", call. = FALSE)
    }

    table <- df[, c(decoy, score)]
    names(table) <- c("decoy", "score")
    table <- stats::na.exclude(table)
    table$score <- as.double(table$score)

    # if variable 'score' is a character, change to continuous
    if (is.character(table$score)) {
        table$score <- as.numeric(as.character(table$score))
    }

    # check whether the selected variables are of the correct class
    if (!is.numeric(table$score)) stop("`score` is not numeric.", call. = FALSE)
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


# Shiny app ---------------------------------------------------------------

#' @import shiny
.select <- function(object, decoy = NULL, score = NULL, log = TRUE, nBins = 50) {
    out <- list(selDecoy = decoy, selScore = score, log = log, nBins = nBins)
    fv <- colnames(object)
    on.exit(return(out))
    ui <- fluidPage(
        title = "Evaluate Decoys",
        sidebarLayout(
            sidebarPanel(
                actionButton("stop", "Stop app"),
                checkboxInput("log", "-log10 transform variable?", out$log),
                selectInput("decoyVar", "Select Decoy",
                    as.list(fv),
                    selected = out$selDecoy
                ),
                selectInput("scoreVar", "Select Score",
                    as.list(fv),
                    selected = out$selScore
                ),
                numericInput("nBins", "Number of bins in histogram", value = out$nBins, min = 2, max = 1000)
            ),
            mainPanel(
                tabsetPanel(
                    type = "tabs",
                    tabPanel(
                        "Info", HTML(
                            paste(
                                h4("Use the left panel to"),
                                "<ul>
                        <li>Indicate if the scores should be log-transformed;</li>
                        <li>Select a <i>Decoy</i> variable</li>
                        <li>Select a <i>Score</i> variable</li>
                        <li>Select the number of <i>Bins</i> for the histogram </li>
                        </ul>",
                                h4("The following conditions must be met"),
                                "<ul>
                        <li> Decoy variable is a boolean that indicates if the score belongs to a target or a decoy</li>
                        <li> Score variable contains the scores of the search engine, which have to be continuous (larger scores are assumed to be better. E-values are typically -log10(e-value) transformed.)</li>
                        <li> You can verify the selected data and settings using the table below, and, the plots in the histogram and PP-Plot tabs.</li>
                        </ul>",
                                "Press button <b>Stop app</b> to return to R, and, to generate and store the plot in your R session."
                            )
                        ),
                        h4("Data table with selected variables"),
                        dataTableOutput("fd")
                    ),
                    tabPanel(
                        "Histogram", plotOutput("histPlot",
                            dblclick = "hist_dblclick",
                            brush = brushOpts(id = "hist_brush", resetOnNew = TRUE)
                        ),
                        p("Brush and double-click on a selected area to zoom in on x-coordinate.
                               Double click outside a selected area to zoom out")
                    ),
                    tabPanel(
                        "PP-Plot", plotOutput("PPplot",
                            dblclick = "PPplot_dblclick",
                            brush = brushOpts(id = "PPplot_brush", resetOnNew = TRUE)
                        ),
                        p("Brush and double-click on a selected area to zoom in on y-coordinate.
                               Double click outside a selected area to zoom out")
                    )
                )
            )
        )
    )

    server <- function(input, output) {
        # Terminate app
        observeEvent(input$stop, {
            stopApp(returnValue = out)
        })

        # Make Data Table
        output$fd <- renderDataTable({
            out$selDecoy <<- input$decoyVar
            out$selScore <<- input$scoreVar
            out$log <<- input$log
            out$nBins <<- input$nBins
            object[, c(out$selDecoy, out$selScore), drop = FALSE]
        })

        # Select Data
        data <- reactive({
            data <- object[, c(input$decoyVar, input$scoreVar)]
            data <- na.exclude(data)
            names(data) <- c("decoy", "score")
            data$score <- as.double(data$score)
            if (input$log) data$score <- -log10(data$score)
            data
        })


        ############
        # Histogram#
        ############

        output$histPlot <- renderPlot({
            binwidth <- diff(range(data()$score, na.rm = TRUE)) / input$nBins
            basePlot <- ggplot(
                data(),
                aes(score, fill = decoy, col = I("black"))
            )
            if (is.logical(data()$decoy) & is.numeric(data()$score)) {
                  basePlot +
                      geom_histogram(alpha = 0.5, bins = input$nBins, position = "identity") +
                      labs(x = "Score", y = "Counts", title = paste0("Histogram of targets and decoys\n")) +
                      coord_cartesian(xlim = histRanges$x, ylim = histRanges$y, expand = TRUE) +
                      theme_bw() +
                      theme(
                          plot.title = element_text(size = rel(1.5)),
                          axis.title = element_text(size = rel(1.2)),
                          axis.text = element_text(size = rel(1.2)),
                          axis.title.y = element_text(angle = 0)
                      )
              }
        })

        histRanges <- reactiveValues(x = NULL, y = NULL)

        observeEvent(input$hist_dblclick, {
            brush <- input$hist_brush
            if (!is.null(brush)) {
                  histRanges$x <- c(brush$xmin, brush$xmax)
              } else {
                  histRanges$x <- NULL
              }
        })

        #########

        #########
        # PP-plot#
        #########
        output$PPplot <- renderPlot({
            basePlot <- ggplot()
            if (is.logical(data()$decoy) & is.numeric(data()$score)) {
                  basePlot <- .ppPlot(data(), PPplotRanges$y)[[1]]
              }
            basePlot
        })
        PPplotRanges <- reactiveValues(x = c(0, 1), y = c(0, 1))

        observeEvent(input$PPplot_dblclick, {
            brush <- input$PPplot_brush
            if (!is.null(brush)) {
                  PPplotRanges$y <- c(0, brush$ymax)
              } else {
                  PPplotRanges$y <- c(0, 1)
              }
        })
        ##########
    }
    app <- list(ui = ui, server = server)
    runApp(app)
}
