# Use one of the example files in the mzID package
exampleFile <- system.file('extdata', '55merge_tandem.mzid', package = 'mzID')
mzIDexample <- mzID::mzID(exampleFile)

test_that("decoyScoreTable() works", {
    out <- decoyScoreTable(mzIDexample,
        decoy = "isdecoy", score = "x\\!tandem:expect", log10 = TRUE
    )
    expect_s3_class(out, "data.frame")
    expect_named(out, c("decoy", "score"))
    expect_type(out$decoy, "logical")
    expect_type(out$score, "double")
})

file <- system.file("mzid", "Tandem.mzid.gz", package = "msdata")
mzRexample <- mzR::openIDfile(file)

test_that("decoyScoreTable() works for mzRident objects", {
    out <- decoyScoreTable(mzRexample,
        decoy = "isDecoy", score = "X.Tandem.expect", log10 = TRUE
    )
    expect_s3_class(out, "data.frame")
    expect_named(out, c("decoy", "score"))
    expect_type(out$decoy, "logical")
    expect_type(out$score, "double")

})

df <- data.frame(
    isdecoy = sample(c(TRUE, FALSE), 100, replace = TRUE),
    score = exp(rnorm(100))
)
test_that("decoyScoreTable() works for data.frame objects", {
    out <- decoyScoreTable(df, decoy = "isdecoy", score = "score", log10 = TRUE)
    expect_s3_class(out, "data.frame")
    expect_named(out, c("decoy", "score"))
    expect_type(out$decoy, "logical")
    expect_type(out$score, "double")

    ## Should coerce character 'score' to double
    df$score <- as.character(df$score)
    out <- decoyScoreTable(df, decoy = "isdecoy", score = "score", log10 = TRUE)
    expect_type(out$score, "double")
})

test_that("decoyScoreTable() breaks for wrong variables", {
    expect_error(decoyScoreTable(mzIDexample,
        decoy = "nope", score = "x\\!tandem:expect", log10 = TRUE
    ))
    expect_error(decoyScoreTable(mzIDexample,
        decoy = "isdecoy", score = "nope", log10 = TRUE
    ))
})

test_that("decoyScoreTable() fails for wrong data format", {
    expect_error(
        ## Deliberately used wrong 'decoy' argument
        decoyScoreTable(df, decoy = "score", score = "score", log10 = TRUE)
    )

    ## Object of wrong class
    x <- as.list(df)
    expect_error(
        decoyScoreTable(x, decoy = "isdecoy", score = "score", log10 = TRUE)
    )
})


# evalTargetDecoys() ------------------------------------------------------

test_that("evalTargetDecoys() works", {
    decoyPlots <- evalTargetDecoys(mzIDexample,
        decoy = "isdecoy", score = "x\\!tandem:expect", log10 = TRUE
    )
    expect_type(decoyPlots, "list")
    expect_length(decoyPlots, 5)
})
