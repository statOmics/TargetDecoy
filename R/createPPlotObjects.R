##########################################################
# FUNCTION 2a
##########################################################
### Function to make a list of the processed mzIDFiles ###
##########################################################
# returns a list of tables with decoys (TRUE/FALSE) and scores
# mzObjects: list of read mzID files
# decoy: string,  name of the variable isdecoy in this mzID file (typical "isdecoy" or "isDecoy")
# score: string, name of the variabele score in this mzID file
# log10: boolean, TRUE if the scores from the mzID are returned log10 transformed

processObjects <- function(mzObjects, decoy = NULL, score = NULL, log10 = TRUE) {
    if (!is.null(decoy) & length(decoy) == 1) decoy <- rep(decoy, length(mzObjects))
    if (!is.null(score) & length(score) == 1) score <- rep(score, length(mzObjects))
    if (!is.null(log10) & length(log10) == 1) log10 <- rep(log10, length(mzObjects))

    decoyScoreTables <- lapply(
        1:length(mzObjects),
        function(i, decoy, score, log10, mzObjects) {
              decoyScoreTable(mzObjects[[i]], decoy[i], score[i], log10[i])
          },
        decoy = decoy,
        score = score,
        log10 = log10,
        mzObjects = mzObjects
    )
    names(decoyScoreTables) <- names(mzObjects)
    return(decoyScoreTables)
}


#################################################################################
# FUNCTION 4
#################################################################################
############## Function for ggplot2-like colour scale in HCL space ##############
#################################################################################

gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


#########################################################################################
#########################################################################################
# Function to create all the PP plots in one figure for scores from files from one engine
#########################################################################################
#########################################################################################
# mzObjects: list of mzID files
# decoy: string,  name of the variable isdecoy in this mzID file (typical "isdecoy" or "isDecoy")
# score: string, name of the variabele score in this mzID file
# log10: boolean, TRUE if the scores from the mzID are returned log10 transformed
# sameSearchEngine: boolean, TRUE if scores result from the same search engine

createPPlotObjects <- function(mzObjects, decoy = NULL, score = NULL, log10 = TRUE, sameSearchEngine = TRUE) {
    if (!class(mzObjects) == "list") mzObjects <- list(mzObjects)

    if ((is.null(decoy) | is.null(score) | is.null(log10)) & sameSearchEngine) {
        object <- mzObjects[[1]]
        df <- .getDF(object)
        if (missing(log10)) log10 <- TRUE
        out <- .select(df, decoy, score, log10)
        decoy <- out$selDecoy # categorical
        score <- out$selScore # continu
        log10 <- out$log
    }

    decoyScoreTables <- processObjects(mzObjects, decoy, score, log10)

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
    ) %>% bind_rows()

    # Plot with original pi0 line and points
    plot1 <- ppData %>%
        ggplot(aes(Fdp, Ftp, color = search)) +
        geom_point() +
        ggtitle("PP plot") +
        ylab("FTarget") +
        xlab("FDecoy") +
        theme_bw() +
        geom_abline(
            slope = ppData %>% pull(pi0) %>% unique(),
            color = gg_color_hue(length(mzObjects))
        )

    # Plot with rescaled pi0 line on 0 and points
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
