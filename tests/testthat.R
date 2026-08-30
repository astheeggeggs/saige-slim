# Standard testthat runner. Under plain R CMD check this runs only the fast
# tier (no SAIGE_TEST_LEVEL set -> helper defaults to "fast" when invoked via
# check). The intended entry point is harness/run_tests.sh, which pins the
# determinism env and sets SAIGE_TEST_LEVEL.
library(testthat)
library(SAIGE)
test_check("SAIGE")
