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


# Function for ggplot2-like colour scale in HCL space
gg_color_hue <- function(n) {
    hues <- seq(15, 375, length.out = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}
