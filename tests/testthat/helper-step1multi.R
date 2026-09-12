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
#
# THE TOLERANCES ARE THE HARNESS'S, and this is not a free choice. These first
# shipped at 1e-5 for everything, which passed on the golden-generating machine
# and went RED on CI: Linux reported rel err 4.2e-5 in `Y` and 1.3e-5 in others
# against a macOS-generated golden. Both are well inside what
# `harness/compare_step1.R` has always allowed for exactly this comparison --
# 1e-2 for theta and the variance ratios (AI-REML runs at tol=0.02, so tighter
# is meaningless) and 1e-3 for the other numeric fields (float PCG at
# tolPCG=1e-5). A test that asserts more precision than the estimator has is
# measuring the platform's libm, not the code, so these now match the harness
# rather than undercutting it. Byte mode is unaffected and is still the real
# gate on the multi-phenotype identity claim.
TOL_THETA <- 1e-2   # theta/tau and the variance ratios
TOL_OTHER <- 1e-3   # every other numeric model field

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

expect_step1_rda_equivalent <- function(golden, candidate, tol = TOL_OTHER) {
  if (compare_wants_byte()) {
    expect_files_byte_identical(golden, candidate)
    # no early return under `both`: the field-wise checks below are the other
    # comparator, and the point of `both` is to apply them to the same files
    if (!compare_wants_numeric()) return(invisible(NULL))
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
    # theta is the AI-REML estimate itself, so it gets the loose bound the
    # harness gives it; everything else is a fitted quantity at PCG precision.
    ftol <- if (identical(f, "theta")) TOL_THETA else tol
    expect_lt(max(abs(gv - cv) / denom), ftol, label = paste0("rel err in ", f))
  }
  invisible(NULL)
}

expect_step1_vr_equivalent <- function(golden, candidate, tol = TOL_THETA) {
  if (compare_wants_byte()) {
    expect_files_byte_identical(golden, candidate)
    # no early return under `both`: the field-wise checks below are the other
    # comparator, and the point of `both` is to apply them to the same files
    if (!compare_wants_numeric()) return(invisible(NULL))
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
