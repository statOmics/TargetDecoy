
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TargetDecoy

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check-bioc](https://github.com/statOmics/TargetDecoy/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/statOmics/TargetDecoy/actions)
[![Codecov test
coverage](https://codecov.io/gh/statOmics/TargetDecoy/branch/master/graph/badge.svg)](https://codecov.io/gh/statOmics/TargetDecoy?branch=master)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/TargetDecoy.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/TargetDecoy)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/TargetDecoy.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/TargetDecoy)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/TargetDecoy.svg)](http://bioconductor.org/packages/stats/bioc/TargetDecoy/)
[![Bioc
support](https://bioconductor.org/shields/posts/TargetDecoy.svg)](https://support.bioconductor.org/tag/TargetDecoy)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/TargetDecoy.svg)](https://bioconductor.org/packages/release/bioc/html/TargetDecoy.html#since)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/TargetDecoy.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/TargetDecoy/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/TargetDecoy.svg)](https://bioconductor.org/packages/release/bioc/html/TargetDecoy.html#since)
<!-- badges: end -->

The goal of **TargetDecoy** is to to generate diagnostic plots to
evaluate the quality of the target decoy approach (TDA).

A first step in the data analysis of Mass Spectrometry (MS) based
proteomics data is to identify peptides and proteins. With this respect
the huge number of experimental mass spectra typically have to be
assigned to theoretical peptides derived from a sequence database.
Search engines are used for this purpose. These tools compare each of
the observed spectra to all candidate theoretical spectra derived from
the sequence data base and calculate a score for each comparison. The
observed spectrum is then assigned to the theoretical peptide with the
best score, which is also referred to as the peptide to spectrum match
(PSM). It is of course crucial for the downstream analysis to evaluate
the quality of these matches. Therefore False Discovery Rate (FDR)
control is used to return a reliable list PSMs. The FDR, however,
requires a good characterisation of the score distribution of PSMs that
are matched to the wrong peptide (bad target hits). In proteomics, the
target decoy approach (TDA) is typically used for this purpose. The TDA
method matches the spectra to a database of real (targets) and nonsense
peptides (decoys). A popular approach to generate these decoys is to
reverse the target database. Hence, all the PSMs that match to a decoy
are known to be bad hits and the distribution of their scores are used
to estimate the distribution of the bad scoring target PSMs. A crucial
assumption of the TDA is that the decoy PSM hits have similar properties
as bad target hits so that the decoy PSM scores are a good simulation of
the target PSM scores. Users, however, typically do not evaluate these
assumptions. To this end we developed TargetDecoy to generate diagnostic
plots to evaluate the quality of the target decoy method.

## Installation

You can install
*[TargetDecoy](https://bioconductor.org/packages/3.16/TargetDecoy)* from
[*Bioconductor*](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("TargetDecoy")
```

The latest development version of **TargetDecoy** can also be installed
from [GitHub](https://github.com/statOmics/TargetDecoy) with:

``` r
BiocManager::install("statOmics/TargetDecoy")
```

## Getting started

Check the
[vignette](https://bioconductor.org/packages/release/bioc/vignettes/TargetDecoy/inst/doc/TargetDecoy.html).
