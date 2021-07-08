#' Select variables
#'
#' @title Select variables of interest
#'
#' @param object A `dataframe`, `mzID` or `mzR` file.
#'
#' @param decoy `boolean` to indicate if the peptide matches to a target or to a decoy.
#'
#' @param score `numeric` indicating the score of the peptide match, obtained by the search engine.
#'
#' @param log10 `boolean` to indicate if the score should be -log10-transformed.
#'
#' @param nBins `numeric` indicating the number of bins in the histogram.
#'
#' @return Diagnostic plots to check the TDA assumptions.
#'
#' @author

.select <- function(object, decoy = NULL, score = NULL, log = TRUE, nBins = 50) {
    require(shiny)
    require(ggplot2)
    out <- list(selDecoy = decoy, selScore = score, log = log, nBins = nBins)
    fv <- colnames(object)
    on.exit(return(out))
    ui <- shiny::fluidPage(
        title = "Evaluate Decoys",
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::actionButton("stop", "Stop app"),
                checkboxInput("log", "-log10 transform variable?", out$log),
                shiny::selectInput("decoyVar", "Select Decoy",
                    as.list(fv),
                    selected = out$selDecoy
                ),
                shiny::selectInput("scoreVar", "Select Score",
                    as.list(fv),
                    selected = out$selScore
                ),
                shiny::numericInput("nBins", "Number of bins in histogram", value = out$nBins, min = 2, max = 1000)
            ),
            shiny::mainPanel(
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
                        shiny::dataTableOutput("fd")
                    ),
                    tabPanel(
                        "Histogram", shiny::plotOutput("histPlot",
                            dblclick = "hist_dblclick",
                            brush = brushOpts(id = "hist_brush", resetOnNew = TRUE)
                        ),
                        p("Brush and double-click on a selected area to zoom in on x-coordinate.
                               Double click outside a selected area to zoom out")
                    ),
                    tabPanel(
                        "PP-Plot", shiny::plotOutput("PPplot",
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
        shiny::observeEvent(input$stop, {
            shiny::stopApp(returnValue = out)
        })

        # Make Data Table
        output$fd <- shiny::renderDataTable({
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
    shiny::runApp(app)
}
