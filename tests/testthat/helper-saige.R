# Shared helpers for the SAIGE test suite. testthat sources helper-*.R
# automatically before running test files.
#
# Tiers (SAIGE_TEST_LEVEL): fast (default) < standard < full — see
# harness/run_tests.sh. Comparison mode (SAIGE_COMPARE_MODE): "byte" (default;
# valid only on the platform that generated harness/goldens) or "numeric"
# (parsed tables, relative tolerance; for cross-platform CI).

suppressPackageStartupMessages({
  library(SAIGE)
  library(data.table)
})

saige_repo <- local({
  p <- normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
  if (!file.exists(file.path(p, "DESCRIPTION"))) p <- Sys.getenv("SAIGE_REPO", unset = p)
  p
})
saige_extdata <- file.path(saige_repo, "extdata")
saige_goldens <- file.path(saige_repo, "harness", "goldens")

saige_level <- function() {
  lv <- Sys.getenv("SAIGE_TEST_LEVEL", unset = "fast")
  if (!lv %in% c("fast", "standard", "full")) stop("bad SAIGE_TEST_LEVEL: ", lv)
  lv
}

skip_unless_level <- function(need) {
  rank <- c(fast = 1, standard = 2, full = 3)
  if (rank[[saige_level()]] < rank[[need]])
    testthat::skip(paste0("requires SAIGE_TEST_LEVEL=", need, " (current: ", saige_level(), ")"))
}

skip_if_no_goldens <- function() {
  testthat::skip_if_not(dir.exists(saige_goldens), "harness/goldens not present")
}

# byte | numeric | both.
#
# `both` exists because the two modes were never testing two different
# GENERATIONS -- they ran the identical step-2 CLI invocations twice and
# differed only in how the outputs were COMPARED. Generation is what costs
# (7 minutes a mode); comparison is reading small tables. So `both` generates
# once and applies both comparators, which is what the gate wanted all along.
#
# It is not a substitute for CI. Byte-exactness is only a valid claim on the
# platform that generated the goldens, so CI stays on `numeric` with its OWN
# generation -- the cross-platform check is a different question from this one.
#
# VALIDATED, because the old reader silently treated every value except "byte"
# as numeric: SAIGE_COMPARE_MODE=bytes ran the weaker comparator and said
# nothing. A typo must not quietly reduce coverage.
saige_compare_mode <- function() {
  m <- Sys.getenv("SAIGE_COMPARE_MODE", unset = "byte")
  if (!m %in% c("byte", "numeric", "both")) {
    stop("bad SAIGE_COMPARE_MODE: '", m, "' (byte, numeric or both)")
  }
  m
}
# Which assertions the current mode licenses. A comparison site asks these
# rather than testing the mode string, so `both` cannot be forgotten at one
# site and silently degrade to numeric there.
compare_wants_byte <- function() saige_compare_mode() %in% c("byte", "both")
compare_wants_numeric <- function() saige_compare_mode() %in% c("numeric", "both")

# Symmetric relative error, same formula as harness/compare_step1.R.
#
# Non-finite values are compared exactly, not by the ratio: `Inf - Inf` is NaN,
# so without this a column containing an infinity could never compare equal to
# itself. Two goldens print them -- plink_marker_sparse.txt and
# plink_marker_noadjCov.txt carry BETA = -inf / SE = nan on the two fixture
# markers that are collinear with the covariates (LEDGER #44) -- so this is a
# live path, not a hypothetical. A non-finite value facing a finite one, or an
# infinity of the opposite sign, is a maximal error rather than a match.
rel_err <- function(a, b, tol = 1e-5) {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b)) return(Inf)
  d <- abs(a - b) / (abs(a) + abs(b) + .Machine$double.eps)
  d[is.na(a) & is.na(b)] <- 0                       # NA and NaN both land here
  same_inf <- is.infinite(a) & is.infinite(b) & (sign(a) == sign(b))
  d[same_inf] <- 0
  if (any(is.finite(a) != is.finite(b))) return(Inf)  # finite vs not
  if (anyNA(d)) return(Inf)  # NA in one but not the other

  # ABSOLUTE FLOOR, RELATIVE TO THE COLUMN'S OWN SCALE (LEDGER #44, second half).
  # A relative test cannot compare two values that are both essentially zero:
  # `Tstat` on a covariate-collinear marker is ~3e-15 where the same column runs
  # to 14.5, and the two platforms' summation orders give different noise, so
  # d ~ 1 for a pair that carries no information whatsoever. Treat a pair as
  # matching when BOTH sides are below tol * max|column|, i.e. when the values
  # are more than five orders of magnitude below the largest value they sit
  # beside. This cannot hide a real difference: anything large enough to affect
  # a conclusion is, by construction, above the floor and still compared
  # relatively. It is scale-free -- a column of genuinely tiny numbers has a
  # correspondingly tiny floor -- and it is NOT a per-row exemption list.
  scale <- suppressWarnings(max(abs(c(a, b))[is.finite(c(a, b))], na.rm = TRUE))
  if (is.finite(scale) && scale > 0) {
    negligible <- is.finite(a) & is.finite(b) &
      pmax(abs(a), abs(b)) <= tol * scale
    d[negligible] <- 0
  }
  max(d, 0)
}

# Run one of the extdata CLI wrappers (step1_fitNULLGLMM.R / step2_SPAtests.R /
# createSparseGRM.R) with cwd = extdata/ using the same R as this session
# (which run_tests.sh launched under pixi with the determinism env pinned).
# Returns the exit status; stdout/stderr go to `log`.
# `env` sets variables for the CHILD only, as name = value pairs. Needed because
# a few behaviours are only reachable through the environment -- notably
# SAIGE_PLINK_NO_MMAP, which forces the .bed pread fallback (S5). Set via
# system2()'s own `env`, so nothing leaks into this process or into other
# parallel test workers.
run_saige_cli <- function(script, args, log = tempfile(fileext = ".log"),
                          env = character()) {
  rscript <- file.path(R.home("bin"), "Rscript")
  owd <- setwd(saige_extdata); on.exit(setwd(owd), add = TRUE)
  envv <- if (length(env)) paste0(names(env), "=", unname(env)) else character()
  status <- suppressWarnings(system2(rscript, c(script, args), stdout = log,
                                     stderr = log, env = envv))
  attr(status, "log") <- log
  status
}

expect_cli_ok <- function(status) {
  if (!identical(as.integer(status), 0L)) {
    lg <- attr(status, "log")
    tail_txt <- if (!is.null(lg) && file.exists(lg))
      paste(utils::tail(readLines(lg, warn = FALSE), 15), collapse = "\n") else "<no log>"
    testthat::fail(paste0("CLI exited with status ", as.integer(status), "; log tail:\n", tail_txt))
  }
  testthat::succeed()
}

# ---- comparators -----------------------------------------------------------
# Pure predicate layer (no testthat expectations) + thin expect_* wrappers.
# The pure layer exists so test-meta-comparators.R can verify detection power
# directly on return values.

files_byte_identical <- function(golden, candidate) {
  if (!file.exists(golden) || !file.exists(candidate)) return(FALSE)
  identical(unname(tools::md5sum(golden)), unname(tools::md5sum(candidate)))
}

read_saige_table <- function(path) {
  # na.strings = character(0): keep literal "NA" as the string "NA" so the
  # numeric-column classifier and NA-pattern comparison see what was printed.
  utils::read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    check.names = FALSE, colClasses = "character", quote = "",
                    comment.char = "", na.strings = character(0))
}

# Parse a printed column, including the non-finite values.
#
# The spelling of a non-finite is the C library's, not R's: glibc prints a NaN
# whose sign bit is set as "-nan" where macOS libc prints "nan", for the same
# computation. Goldens P27 and P30 contain them (rs20/rs22, whose covariate-
# adjusted variance underflows -- LEDGER #44), so this is a live path.
#
# as.numeric() handles "inf"/"-inf" but returns NaN for "nan"/"-nan", which
# anyNA() cannot tell apart from a parse FAILURE. That is why a column holding a
# NaN used to be classified non-numeric and compared as text -- and text is
# exactly what differs across platforms. Map the spellings explicitly.
.num <- function(x) {
  x <- trimws(x)
  out <- suppressWarnings(as.numeric(x))
  out[grepl("^[+-]?nan$", x, ignore.case = TRUE)] <- NaN
  out
}

# Column is treated as numeric if every non-"NA" entry parses as a number or as
# a non-finite. A real parse failure gives NA (not NaN), which is what fails.
.is_numeric_col <- function(x) {
  x2 <- x[!(x %in% c("NA", ""))]
  if (length(x2) == 0) return(TRUE)
  v <- .num(x2)
  all(!is.na(v) | is.nan(v))
}

# Tolerance-based table comparison for step-2 outputs: structure and string
# columns must match exactly; numeric columns within `tol` symmetric relative
# error ("NA" must match "NA"). Returns list(ok, problems).
#
# `drop_markers` removes named rows from BOTH tables after the structural checks
# (column sets, dims and the MarkerID column itself are still compared over the
# full tables, so a run that lost or reordered a marker still fails). It exists
# for rows whose statistics are numerically meaningless and therefore not
# portable across platforms -- see the LEDGER #44 note at its only call site.
# Use it for that, and never to quiet a difference that carries information.
tables_equivalence <- function(golden, candidate, tol = 1e-5,
                               drop_markers = character()) {
  problems <- character()
  if (!file.exists(golden)) problems <- c(problems, paste("missing golden:", golden))
  if (!file.exists(candidate)) problems <- c(problems, paste("missing candidate:", candidate))
  if (length(problems)) return(list(ok = FALSE, problems = problems))
  g <- read_saige_table(golden); c_ <- read_saige_table(candidate)
  if (!identical(colnames(g), colnames(c_)))
    return(list(ok = FALSE, problems = paste("column sets differ:", basename(golden))))
  if (!identical(dim(g), dim(c_)))
    return(list(ok = FALSE, problems = paste("dims differ:", basename(golden))))
  if (length(drop_markers)) {
    if (is.null(g[["MarkerID"]]))
      return(list(ok = FALSE, problems = paste("drop_markers needs a MarkerID column:",
                                               basename(golden))))
    if (!identical(g[["MarkerID"]], c_[["MarkerID"]]))
      return(list(ok = FALSE, problems = paste("string column differs [MarkerID]:",
                                               basename(golden))))
    if (!all(drop_markers %in% g[["MarkerID"]]))
      return(list(ok = FALSE, problems = paste("drop_markers absent from table:",
                                               basename(golden))))
    keep <- !(g[["MarkerID"]] %in% drop_markers)
    g <- g[keep, , drop = FALSE]
    c_ <- c_[keep, , drop = FALSE]
  }
  for (col in colnames(g)) {
    if (.is_numeric_col(g[[col]]) && .is_numeric_col(c_[[col]])) {
      if (!identical(g[[col]] %in% c("NA", ""), c_[[col]] %in% c("NA", "")))
        problems <- c(problems, paste0("NA pattern differs [", col, "]"))
      err <- rel_err(.num(g[[col]]), .num(c_[[col]]), tol = tol)
      if (!is.finite(err) || err > tol)
        problems <- c(problems, sprintf("numeric [%s] rel_err=%.3e > tol=%.0e", col, err, tol))
    } else if (!identical(g[[col]], c_[[col]])) {
      problems <- c(problems, paste0("string column differs [", col, "]"))
    }
  }
  list(ok = length(problems) == 0, problems = problems)
}

expect_files_byte_identical <- function(golden, candidate) {
  testthat::expect_true(files_byte_identical(golden, candidate),
                        label = paste0("byte-identical: ", basename(golden),
                                       " (exists: g=", file.exists(golden),
                                       " c=", file.exists(candidate), ")"))
}

expect_tables_equivalent <- function(golden, candidate, tol = 1e-5,
                                     drop_markers = character()) {
  r <- tables_equivalence(golden, candidate, tol = tol, drop_markers = drop_markers)
  testthat::expect_true(r$ok, label = paste0("tables equivalent: ", basename(golden),
                                             if (!r$ok) paste0(" [", paste(r$problems, collapse = "; "), "]")))
}

# Mode-dispatching comparison used by the end-to-end tests. `drop_markers`
# reaches only the numeric path: byte mode is a claim about one code path re-run
# on ONE platform, where even a meaningless row reproduces exactly, so nothing is
# excluded there.
expect_output_matches_golden <- function(golden, candidate, tol = 1e-5,
                                         drop_markers = character()) {
  if (compare_wants_byte()) {
    expect_files_byte_identical(golden, candidate)
  }
  if (compare_wants_numeric()) {
    expect_tables_equivalent(golden, candidate, tol = tol,
                             drop_markers = drop_markers)
  }
}
