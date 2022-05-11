#' @import shiny
#' @import miniUI
.select_vars <- function(object, decoy = NULL, score = NULL,
                         log = TRUE, nBins = 50) {

    df <- .getDF(object)

    score_choices <- .find_vars(df, is.numeric)
    decoy_choices <- .find_vars(df, is.logical)

    ui <- miniPage(
        gadgetTitleBar("TargetDecoy Gadget"),
        miniTabstripPanel(
            miniTabPanel("Variables", icon = icon("sliders-h"),
                miniContentPanel(
                    checkboxInput("log", "-log10 transform variable?",
                        value = log
                    ),
                    selectInput("decoyVar", "Select Decoy",
                        choices = decoy_choices,
                        selected = decoy
                    ),
                    selectInput("scoreVar", "Select Score",
                        choices = score_choices,
                        selected = score
                    ),
                    numericInput("nBins", "Number of bins in histogram",
                        value = nBins, min = 2, max = 1000
                    )
                )
            ),
            miniTabPanel("Histogram",
                miniContentPanel(
                    plotOutput("hist", height = "100%")
                )

            ),
            miniTabPanel("PP-plot",
                miniContentPanel(
                    plotOutput("PPplot", height = "100%")
                )
            ),
            miniTabPanel("Data", icon = icon("table"),
                miniContentPanel(
                    dataTableOutput("data")
                )
            )
        )
    )

    server <- function(input, output, session) {
        vars <- reactive({
            validate(
                need(is.logical(df[[input$decoyVar]]),
                    "`decoy` variable should be logical"),
                need(is.numeric(df[[input$scoreVar]]),
                    "`score` variable should be numeric.")
            )

            list(
                decoy = input$decoyVar,
                score = input$scoreVar,
                log = input$log,
                nBins = input$nBins
            )
        })

        observeEvent(input$done, {
            stopApp(returnValue = vars())
        })

        output$hist <- renderPlot({
            vars <- vars()
            evalTargetDecoysHist(df,
                decoy = vars$decoy, score = vars$score,
                log10 = vars$log, nBins = vars$nBins
            )
        })

        ## FIXME: validate that at least some decoys are present
        output$PPplot <- renderPlot({
            vars <- vars()
            evalTargetDecoysPPPlot(df,
                decoy = vars$decoy, score = vars$score, log10 = vars$log
            )
        })

        output$data <- renderDataTable({
            vars <- vars()
            decoyScoreTable(df,
                decoy = vars$decoy, score = vars$score, log10 = vars$log
            )
        })
    }

    runGadget(ui, server)
}


## Helper to find variables from a specific type
.find_vars <- function(data, filter) {
    names(data)[vapply(data, filter, logical(1))]
}
