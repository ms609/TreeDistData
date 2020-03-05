## Reminders when releasing for CRAN
release_questions <- function() {
  c(
    "Is the code free of #TODOs?",
    "Have you updated REFERENCES.bib, inst/CITATION and README.md with a citation to the published study?",
    "Have you updated the version number in NEWS & DESCRIPTION?",
    "Have you cleared GitHub issues for this release milestone?",
    "Have you checked the Vignettes for sanity?"
  )
}

# Additional tests:
#
# check_win_devel(); check_rhub()
# #rhub::check_on_windows()
# rhub::check_with_rdevel() # redundifies check_on_debian()
# rhub::check_on_ubuntu()
# rhub::check_on_fedora()
# rhub::check_on_centos()
# rhub::check_with_valgrind() # runs the build and check on Linux, in valgrind to find memory leaks and pointer errors.
# rhub::check_with_sanitizers() # runs all package package tests, examples and vignettes with Address Sanitizer and Undefined Behavior Sanitizer.
# list_my_checks() # list_package_checks
# revdepcheck::revdep_check()
# build_vignettes()
# tools::resaveRdaFiles('data', compress='auto') - is default of bzip2 the optimal?
# tools::checkRdaFiles('data') - set optimal compression in `data-raw`
