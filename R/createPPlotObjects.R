#' Create all the PP plots in one figure for scores from multiple objects
#'
#' @param object_list List of \linkS4class{mzID} or \linkS4class{mzRident}
#'   objects. If named, the names will be used in the legend of the plot. If
#'   not, names will be extracted from the data files in case of *mzID* or
#'   *mzRident* objects.
#' @inheritParams decoyScoreTable
#'
#' @return
#' One PP plot with all original pi0, and a standardized / rescaled PP plot with
#' all `pi0` set to 0.
#'
#' @author Elke Debrie, Lieven Clement
#'
#' @export
#'
#' @examples
#' library(mzID)
#'
#' ## Use two example files from the mzID package
#' exampleFiles <- system.file(
#'     "extdata", c("55merge_omssa.mzid", "55merge_tandem.mzid"),
#'     package = "mzID"
#' )
#' mzObjects <- lapply(exampleFiles, mzID)
#'
#' createPPlotObjects(mzObjects,
#'     decoy = "isdecoy",
#'     score = c("omssa:evalue", "x\\!tandem:expect"),
#'     log10 = TRUE
#' )
createPPlotObjects <- function(object_list,
                               decoy = NULL, score = NULL, log10 = TRUE) {
    if (!is.list(object_list)) object_list <- list(object_list)

    tables <- processObjects(
        object_list = object_list,
        decoy = decoy, score = score, log10 = log10
    )

    ppData <- ppScoresData(tables)
    ppScoresPlots(ppData)
}


# Helper to make decoy score tables from multiple mzID objects
processObjects <- function(object_list, decoy, score, log10) {
    arg_list <- .check_args(
        object_list = object_list,
        decoy = decoy, score = score, log10 = log10
    )
    object_list <- arg_list$object_list
    decoy <- arg_list$decoy
    score <- arg_list$score
    log10 <- arg_list$log10

    out <- vector("list", length = length(object_list))
    for (i in seq_along(object_list)) {
        out[[i]] <- decoyScoreTable(
            object = object_list[[i]],
            decoy = decoy[[i]],
            score = score[[i]],
            log10 = log10[[i]]
        )
    }
    names(out) <- names(object_list)

    ## If no names given, use raw file names from mz objects
    if (is.null(names(out))) {
        names(out) <- vapply(object_list, .get_object_name, character(1))
    }

    out
}


.get_object_name <- function(object) {
    if (is(object, "mzID")) {
        fname <- mzID::files(object)$id
        out <- basename(fname)
    } else if (is(object, "mzRident")) {
        fname <- mzR::fileName(object)
        out <- basename(fname)
    } else {
        out <- NULL
    }
    out
}
