# Use two example files from the mzID package
exampleFiles <- system.file(
    "extdata", c("55merge_omssa.mzid", "55merge_tandem.mzid"),
    package = "mzID"
)
mzObjects <- lapply(exampleFiles, mzID::mzID)

test_that("createPPlotObjects() works", {
    out <- createPPlotObjects(mzObjects,
        decoy = "isdecoy",
        score = c("omssa:evalue", "x\\!tandem:expect"),
        log10 = TRUE
    )

    expect_type(out, "list")
    expect_length(out, 2)
})

test_that("createPPlotObjects() works for single objects", {
    out <- createPPlotObjects(mzObjects[[1]],
        decoy = "isdecoy", score = "omssa:evalue", log10 = TRUE
    )

    expect_type(out, "list")
    expect_length(out, 2)
})

test_that("createPPlotObjects() fails for incompatible arg lengths", {
    expect_error(
        createPPlotObjects(mzObjects,
            score = c("omssa:evalue", "x\\!tandem:expect"),
            decoy = "isdecoy", log10 = c(TRUE, FALSE, TRUE)
        )
    )
    expect_error(
        createPPlotObjects(mzObjects,
            score = c("omssa:evalue", "x\\!tandem:expect"),
            decoy = rep("isdecoy", 3), log10 = TRUE
        )
    )
})

## Add some more examples to the list in different formats
mzRexample <- mzR::openIDfile(system.file("mzid", "Tandem.mzid.gz", package = "msdata"))
mzObjects[[3]] <- mzRexample
mzObjects[[4]] <- data.frame(
    decoy = sample(c(TRUE, FALSE), 100, replace = TRUE), score = exp(rnorm(100))
)

test_that("createPPlotObjects() works for various object formats", {
    out <- createPPlotObjects(mzObjects,
        decoy = c("isdecoy", "isdecoy", "isDecoy", "decoy"),
        score = c("omssa:evalue", "x\\!tandem:expect", "X.Tandem.expect", "score"),
        log10 = TRUE
    )

    expect_type(out, "list")
    expect_length(out, 2)
})
