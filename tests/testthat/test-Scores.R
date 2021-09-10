## Use one of the example files in the mzID package
exampleFile <- system.file("extdata", "55merge_tandem.mzid", package = "mzID")
mzIDexample <- mzID::mzID(exampleFile)

test_that("createPPlotScores() works", {
    out <- createPPlotScores(mzIDexample,
        scores = c("x\\!tandem:hyperscore", "x\\!tandem:expect"),
        decoy = "isdecoy", log10 = TRUE
    )

    expect_type(out, "list")
    expect_length(out, 2)
})

test_that("createPPLotScores() fails for incompatible arg lengths", {
    expect_error(
        createPPlotScores(mzIDexample,
            scores = c("x\\!tandem:hyperscore", "x\\!tandem:expect"),
            decoy = "isdecoy", log10 = c(TRUE, FALSE, TRUE)
        ),
        "They should have the same length or length == 1."
    )
})
