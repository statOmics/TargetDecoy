################################################
# FUNCTION 1
################################################
### Function for preprocessing the mzIDfiles ###
################################################

# returns a table with decoys (TRUE/FALSE) and scores from a given object
# object: string, read mzID file
# decoy: string,  name of the variable isdecoy in this mzID file (typical "isdecoy" or "isDecoy")
# score: string, name of the variabele score in this mzID file
# log10: boolean, TRUE if the scores from the mzID are returned log10 transformed

decoyScoreTable <- function(object, decoy = NULL, score = NULL, log10 = TRUE) {
    if (class(object) == "mzID") {
        dataFrame <- flatten(object)
    } else if (class(object) == "data.frame") {
        dataFrame <- object
    } else if (class(object) == "mzRident") {
        dataFrame <- cbind(psms(object), score(object)[, -1])
    } else {
        "object should be of the class mzID, mzRident or dataframe"
    }


    if (is.null(decoy) | is.null(score) | is.null(log10)) {
        if (missing(log10)) log10 <- TRUE
        out <- .select(dataFrame, decoy, score, log10)
        decoy <- out$selDecoy # categorical
        score <- out$selScore # continu
        log10 <- out$log
    }

    table <- dataFrame[, c(decoy, score)]
    names(table) <- c("decoy", "score")
    table <- na.exclude(table)
    table$score <- as.double(table$score)

    # if variable 'score' is a character, change to continuous
    if (class(table$score) == "character") {
        table$score <- as.numeric(as.character(table$score))
    }

    # perform log10-transformation on variable 'score' if so indicated
    if (log10) {
        table$score <- -log10(as.numeric(table$score))
    }

    return(table)
}

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
# FUNCTION 3
#################################################################################
### Function to make the basic plot that is used to plot one ore more ppPlots ###
#################################################################################
# basic plot to which all other points will be added
# startplot: basic plot to start with (typical an empty plot)
# title: title for the plot
# ylab: label on the y-axis

initializePlotDeco <- function(startplot, title = NULL, ylab = NULL) {
    return(startplot +
        labs(title = title) +
        coord_cartesian(xlim = c(0, 1), ylim = NULL, expand = TRUE) +
        xlab(expression(paste("FDecoy"))) +
        ylab(ylab) +
        theme_bw() +
        theme(
            plot.title = element_text(size = rel(1.5)),
            axis.title = element_text(size = rel(1.2)),
            axis.text = element_text(size = rel(1.2)),
            axis.title.y = element_text(angle = 0)
        ))
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
        if (class(object) == "mzID") {
            dataFrame <- flatten(object)
        } else if (class(object) == "data.frame") {
            dataFrame <- object
        } else if (class(object) == "mzRident") {
            dataFrame <- cbind(psms(object), score(object)[, -1])
        } else {
            "object should be of the class mzID, mzRident or dataframe"
        }
        if (missing(log10)) log10 <- TRUE
        out <- .select(dataFrame, decoy, score, log10)
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
