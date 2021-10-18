#' Prepare score table for decoys
#'
#' Takes an input object and returns a score table for the decoys. If one of the
#' arguments besides `object` is missing, will open a Shiny app to interactively
#' select the variables.
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
decoyScoreTable <- function(object, decoy = NULL, score = NULL, log10 = TRUE) {
    df <- .getDF(object)

    if (is.null(decoy) || is.null(score) || is.null(log10)) {
        if (missing(log10)) log10 <- TRUE
        out <- .select(df, decoy, score, log10)
        decoy <- out$selDecoy # categorical
        score <- out$selScore # continu
        log10 <- out$log
    }
    
    ###############
    # VERVANGEN?
    
    if (!(score %in% colnames(df)) | !(decoy %in% colnames(df))) {
        stop("Column(s) not found:", call. = FALSE)
    }
    
    ###############

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


