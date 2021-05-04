## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Do you plan to add a reference to Smith (2021)?",
    "Have you updated REFERENCES.bib, inst/CITATION and README.md with a citation to the published study?",
    "Have you updated the version number in NEWS & DESCRIPTION?",
    "Have you checked the Vignettes for sanity?"
  )
}

# Additional steps:
#
# Propogate changes in README.md to R/TreeDistData-package.R


# Additional tests:
#
# spell_check()
# pkgdown::build_reference_index()
#
# run_examples()
# build_vignettes()
#
# devtools::check_win_devel(); rhub::check_for_cran()
# rhub::check_with_valgrind() # runs the build and check on Linux, in valgrind to find memory leaks and pointer errors.
# rhub::check_with_sanitizers() # runs all package package tests, examples and vignettes with Address Sanitizer and Undefined Behavior Sanitizer.
#
# revdepcheck::revdep_check()
#
# codemetar::write_codemeta()
#
# tools::resaveRdaFiles('data', compress = 'auto') - is default bzip2 the optimal?
# tools::checkRdaFiles('data') - set optimal compression in `data-raw`
