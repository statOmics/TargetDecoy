# Function for preprocessing mzIDfiles
# returns a table with decoys (TRUE/FALSE) and scores from a given mzID file

decoyScoreTable <- function(object, decoy = NULL, score = NULL, log10 = TRUE) {

    # dataFrame <- flatten(mzID(paste0(homepath,file))) # ! functie maakt nog gebruik van homepath !
    # dataFrame <- flatten(mzID(file))
    if (is(object, "mzID")) {
        dataFrame <- mzID::flatten(object)
    } else if (is(object, "data.frame")) {
        dataFrame <- object
    } else if (is(object, "mzRident")) {
        dataFrame <- cbind(psms(object), score(object)[, -1])
    } else {
        "object should be of the class mzID, mzRident or dataframe"
    }


    if (is.null(decoy) | is.null(score) | is.null(log10)) {
        if (missing(log10)) log10 <- TRUE
        out <- .selectMini(dataFrame, decoy, score, log10)
        decoy <- out$selDecoy # categorical
        score <- out$selScore # continu
        log10 <- out$log
    }

    # if("isDecoy" %in% colnames(dataFrame)){
    #   decoy <- "isDecoy"
    # } else if ("isdecoy" %in% colnames(dataFrame)){

    # decoy <- "isdecoy"
    # } else
    # {
    #   out <- .selectMiniDecoy(dataFrame,decoy,score,log10)
    #   decoy <- out$selDecoy #categorical
    #   score <- out$selScore #continu
    #   log10 <- out$log
    # }
    #
    # if("x!tandem:expect" %in% colnames(dataFrame)){
    #   score <- "x!tandem:expect"
    # } else  if("ms-gf:specevalue" %in% colnames(dataFrame)){
    #   score <- "ms-gf:specevalue"
    # } else   {
    #   out <- .selectMiniScore(dataFrame,decoy,score,log10)
    #   decoy <- out$selDecoy #categorical
    #   score <- out$selScore #continu
    #   log10 <- out$log
    # }
    table <- dataFrame[, c(decoy, score)]
    names(table) <- c("decoy", "score")
    table <- na.exclude(table)
    table$score <- as.double(table$score)

    # if variable 'score' is a character, change to continuous
    if (is(table$score, "character")) {
        table$score <- as.numeric(as.character(table$score))
    }

    # perform log10-transformation on variable 'score' if so indicated
    if (log10) {
        table$score <- -log10(as.numeric(table$score))
    }

    return(table)
}


# Function to make a list of the processed mzIDFiles
# returns a list of tables with decoys (TRUE/FALSE) and scores from a given mzID file

processScores <- function(mzObject, decoy = NULL, score, log10 = TRUE) {
    if (!is.null(decoy) & length(decoy) == 1) decoy <- rep(decoy, length(score))
    if (!is.null(log10) & length(log10) == 1) log10 <- rep(log10, length(score))

    decoyScoreTables <- lapply(
        1:length(score),
        function(i, decoy, score, log10, mzObjects) {
              decoyScoreTable(mzObject, decoy[i], score[i], log10[i])
          },
        decoy = decoy,
        score = score,
        log10 = log10,
        mzObjects = mzObject
    )
    names(decoyScoreTables) <- score
    return(decoyScoreTables)
}



# Function for ggplot2-like colour scale in HCL space

gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Function to evaluate assumptions of the Target Decoys Approach for multiple search engines
#'
#' @description Function to create diagnostic PP plots in one plot to evaluate the TDA assumptions
#' for multiple search engines. The function provides the possibility to evaluate each of
#' the sub-engines and the overall itself.
#'
#' @param mzObject A `mzID` file processed by mzID().
#'
#' @param decoy `boolean` to indicate if the peptide matches to a target or to a decoy.
#'
#' @param score `numeric` indicating the score of the peptide match, obtained by the search engine.
#'
#' @param log10 `boolean` to indicate if the score should be -log10-transformed.
#'
#' @param sameSearchEngine `boolean` indicating if runs were processed by same search engine.
#'
#' @return One PP plot with all original pi0, and a standardized / rescaled PP plot with all pi0 set to 0.
#'
#' @author Elke Debrie, Lieven Clement
#'
#' @import stats
#' @import ggplot2
#' @import dplyr
#'
#' @importFrom tibble tibble
#'
#' @export
#'
#' @example
#' # Example with 3 search eniges and a combination of the three
#' mzShaker <- mzID(path = my_path)
#' mzObjects <- lapply(mzFiles, mzID) # list of mzID files
#' createPPlotObjects(mzObjects, decoy = "isdecoy", score = "ms-gf:specevalue") # same engine
#' createPPlotObjects(mzObjects, decoy = "isdecoy", score=rep(c("ms-gf:specevalue","x!tandem:expect"),2)) # different engine


createPPlotScores <- function(mzObject, decoy = NULL, scores, log10 = TRUE) {
    decoyScoreTables <- processScores(mzObject, decoy, scores, log10)

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
    ) %>%
        bind_rows() %>%
        arrange(search)

    plot1 <- ppData %>%
        ggplot(aes(Fdp, Ftp, color = search)) +
        geom_point() +
        ggtitle("PP plot") +
        ylab("FTarget") +
        xlab("FDecoy") +
        theme_bw() +
        geom_abline(
            slope = ppData %>% pull(pi0) %>% unique(),
            color = gg_color_hue(length(scores))
        )

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
