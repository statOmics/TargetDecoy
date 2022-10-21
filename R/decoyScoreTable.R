#' Prepare score table for decoys
#'
#' Takes an input object and returns a score table for the decoys.
#'
#' @param object A `data.frame`, \linkS4class{mzID} or \linkS4class{mzRident} object.
#' @param decoy `character`, name of the variable that indicates if the peptide
#'   matches to a target or to a decoy. When this value is missing, a Shiny
#'   gadget is launched to select it interactively.
#' @param score `numeric`, indicating the score of the peptide match, obtained
#'   by the search engine. When this value is missing, a Shiny gadget is
#'   launched to select it interactively.
#' @param log10 `logical` to indicate if the score should be
#'   `-log10`-transformed. Default: `TRUE`. When this value is missing, a Shiny
#'   gadget is launched to select it interactively.
#'
#' @return
#' A `data.frame` with a logical `"decoy"` column and numeric `"scores"`.
#'
#' @author Elke Debrie, Lieven Clement
#' @export
decoyScoreTable <- function(object, decoy, score, log10 = TRUE) {
    stopifnot(is.logical(log10))

    df <- .getDF(object)

    is_missing <- !(c(decoy, score) %in% colnames(df))
    if (any(is_missing)) {
        missing <- paste(c(decoy, score)[is_missing], collapse = ", ")
        stop("Variable(s) not found: ", missing, call. = FALSE)
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
