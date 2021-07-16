#' Create all the PP plots in one figure for scores from multiple objects
#'
#' @param object_list
#' @inheritParams decoyScoreTable
#'
#' @return
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

    # FIXME: find way to set proper names for tables (extract from mz objects?)
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
    out
}
