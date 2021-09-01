## Helper to compute data for PP-plot
#' @importFrom stats ecdf
.ppData <- function(data) {
    pi0 <- sum(data$decoy) / sum(!data$decoy)

    x <- data$score[!data$decoy]
    Ft <- ecdf(x)
    Fd <- ecdf(data$score[data$decoy])
    z <- Ft(x) - pi0 * Fd(x)
    df <- data.frame(Fdp = Fd(x), Ftp = Ft(x), z = z)
    ylimHlp <- mean(Fd(x) == 1)

    list(df = df, pi0 = pi0, ylimHlp = ylimHlp)
}


## Create PP-plot
.ppPlot <- function(data, ylim = NULL) {
    ppData <- .ppData(data)
    df <- ppData$df
    pi0 <- ppData$pi0
    ylimHlp <- ppData$ylimHlp

    ## Avoid R CMD check "no visible binding" warnings
    Fdp <- Ftp <- NULL

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

   list(ppPlot = ppPlot, yLim = c(0, (1 - ylimHlp)))
}


# Function for ggplot2-like colour scale in HCL space
gg_color_hue <- function(n) {
    hues <- seq(15, 375, length.out = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}
