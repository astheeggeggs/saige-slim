# Comparing STEP-1 outputs between two runs.
#
# helper-saige.R's expect_output_matches_golden is a step-2 comparator: in
# numeric mode it parses a TSV with a header. A step-1 `.rda` is a binary
# serialisation and a `.varianceRatio.txt` is a headerless
# `value type category` file, so neither goes through it -- handing an .rda to
# read.table produces "embedded nulls" and "invalid UTF-8", which is what
# happened the first time these tests ran in numeric mode.
#
# Byte mode asserts what these tests actually claim: ONE code path re-run,
# on one machine, produces the same bytes. Numeric mode exists for the
# golden-generating platform's outputs being compared elsewhere, so there the
# `.rda` is compared field by field at a relative tolerance and the variance
# ratio file column by column.

# The step-1 model fields worth comparing: every numeric one that a refit could
# move. Names are checked first, so a field appearing or disappearing fails.
STEP1_NUMERIC_FIELDS <- c("theta", "coefficients", "linear.predictors",
                          "fitted.values", "Y", "residuals", "cov", "y",
                          "X", "offset")

load_modglmm <- function(path) {
  e <- new.env()
  load(path, envir = e)
  get("modglmm", envir = e)
}

expect_step1_rda_equivalent <- function(golden, candidate, tol = 1e-5) {
  if (saige_compare_mode() == "byte") {
    expect_files_byte_identical(golden, candidate)
    return(invisible(NULL))
  }
  g <- load_modglmm(golden)
  c_ <- load_modglmm(candidate)
  expect_identical(sort(names(g)), sort(names(c_)))
  expect_identical(g$traitType, c_$traitType)
  expect_identical(as.character(g$sampleID), as.character(c_$sampleID))
  for (f in intersect(STEP1_NUMERIC_FIELDS, names(g))) {
    gv <- as.numeric(g[[f]])
    cv <- as.numeric(c_[[f]])
    expect_identical(length(gv), length(cv), label = paste0("length of ", f))
    denom <- pmax(abs(gv), abs(cv), .Machine$double.eps)
    expect_lt(max(abs(gv - cv) / denom), tol, label = paste0("rel err in ", f))
  }
  invisible(NULL)
}

expect_step1_vr_equivalent <- function(golden, candidate, tol = 1e-5) {
  if (saige_compare_mode() == "byte") {
    expect_files_byte_identical(golden, candidate)
    return(invisible(NULL))
  }
  rd <- function(p) {
    utils::read.table(p, header = FALSE, stringsAsFactors = FALSE,
                      colClasses = "character")
  }
  g <- rd(golden)
  c_ <- rd(candidate)
  expect_identical(dim(g), dim(c_))
  # Column 1 is the ratio; the rest are the type and category LABELS, which
  # must agree exactly -- a run that put the right number in the wrong bin is
  # not an equivalent run.
  if (ncol(g) > 1) {
    for (j in 2:ncol(g)) expect_identical(g[[j]], c_[[j]])
  }
  gv <- as.numeric(g[[1]])
  cv <- as.numeric(c_[[1]])
  denom <- pmax(abs(gv), abs(cv), .Machine$double.eps)
  expect_lt(max(abs(gv - cv) / denom), tol)
  invisible(NULL)
}

# Dispatch on the extension, so call sites read as "these two step-1 outputs
# agree".
expect_step1_matches <- function(golden, candidate) {
  if (grepl("\\.rda$", golden)) {
    expect_step1_rda_equivalent(golden, candidate)
  } else {
    expect_step1_vr_equivalent(golden, candidate)
  }
}
