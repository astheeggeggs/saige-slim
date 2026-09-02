# Self-tests for the comparison machinery: prove the comparators pass on
# identical input (null control) and FAIL on planted perturbations (positive
# control). Asserts on the pure predicate layer (files_byte_identical /
# tables_equivalence) so detection power is tested directly. If these break,
# every green result from the suite is meaningless.

make_tsv <- function(df) {
  f <- tempfile(fileext = ".txt")
  utils::write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE)
  f
}

base_df <- data.frame(
  CHR = c("1", "1", "1"),
  MarkerID = c("rs1", "rs2", "rs3"),
  BETA = c("-1.09779", "0.512345", "NA"),
  p.value = c("7.940189E-01", "3.100000E-02", "NA"),
  Is.SPA = c("false", "true", "false"),
  stringsAsFactors = FALSE
)

test_that("rel_err behaves", {
  expect_equal(rel_err(c(1, 2), c(1, 2)), 0)
  expect_gt(rel_err(1, 1.1), 0.04)
  expect_equal(rel_err(NA_real_, NA_real_), 0)
  expect_identical(rel_err(1, NA_real_), Inf)
  expect_identical(rel_err(1, c(1, 2)), Inf)
})

test_that("rel_err on non-finite values: matches pass, mismatches are maximal", {
  # null control: a column that is bit-for-bit itself must compare equal even
  # when it holds infinities and NaNs (goldens do -- see LEDGER #44).
  x <- c(-Inf, Inf, NaN, 1.5)
  expect_equal(rel_err(x, x), 0)
  # positive controls: every way of getting it wrong must still be caught
  expect_identical(rel_err(-Inf, Inf), Inf)      # sign flip
  expect_identical(rel_err(-Inf, -1e300), Inf)   # infinity vs a huge finite
  expect_identical(rel_err(NaN, 0), Inf)         # NaN vs a value
  expect_gt(rel_err(c(-Inf, 1.5), c(-Inf, 1.6)), 0.03)  # finite entries still count
})

test_that("byte comparator: null control passes, perturbation fails", {
  a <- make_tsv(base_df)
  expect_true(files_byte_identical(a, make_tsv(base_df)))

  pert <- base_df; pert$p.value[1] <- "7.940188E-01"   # last-digit flip
  expect_false(files_byte_identical(a, make_tsv(pert)))
  expect_false(files_byte_identical(a, "/nope/missing.txt"))
})

test_that("numeric comparator: null control passes, within-tolerance noise passes", {
  a <- make_tsv(base_df)
  expect_true(tables_equivalence(a, make_tsv(base_df), tol = 1e-5)$ok)

  noisy <- base_df; noisy$BETA[1] <- "-1.09779000000012"   # ~1e-13 rel
  expect_true(tables_equivalence(a, make_tsv(noisy), tol = 1e-5)$ok)
})

test_that("numeric comparator: catches every planted defect class", {
  a <- make_tsv(base_df)

  drift <- base_df; drift$BETA[1] <- "-1.09790"            # ~5e-5 rel > tol
  r <- tables_equivalence(a, make_tsv(drift), tol = 1e-5)
  expect_false(r$ok); expect_match(r$problems, "numeric \\[BETA\\]", all = FALSE)

  s <- base_df; s$Is.SPA[1] <- "true"
  expect_false(tables_equivalence(a, make_tsv(s))$ok)

  n <- base_df; n$BETA[3] <- "0.1"                          # NA -> value
  r <- tables_equivalence(a, make_tsv(n))
  expect_false(r$ok); expect_match(r$problems, "NA pattern", all = FALSE)

  inf_df <- base_df; inf_df$BETA[1] <- "-inf"               # the LEDGER #44 shape
  expect_true(tables_equivalence(make_tsv(inf_df), make_tsv(inf_df))$ok)
  flipped <- inf_df; flipped$BETA[1] <- "inf"
  expect_false(tables_equivalence(make_tsv(inf_df), make_tsv(flipped))$ok)

  # The C library's NaN spelling is platform-dependent: glibc writes "-nan"
  # where macOS libc writes "nan" for the same value (this is what turned CI red
  # on P27). The two must compare EQUAL -- a NaN has no meaningful sign -- while
  # a NaN facing a number must not.
  nan_df <- base_df; nan_df$BETA[1] <- "nan"
  neg_nan <- base_df; neg_nan$BETA[1] <- "-nan"
  expect_true(tables_equivalence(make_tsv(nan_df), make_tsv(neg_nan))$ok)
  num <- base_df; num$BETA[1] <- "0.5"
  expect_false(tables_equivalence(make_tsv(nan_df), make_tsv(num))$ok)
  # and a column of genuine non-numbers is still compared as text
  txt <- base_df; txt$BETA[1] <- "banana"
  expect_false(tables_equivalence(make_tsv(nan_df), make_tsv(txt))$ok)

  expect_false(tables_equivalence(a, make_tsv(base_df[, -3]))$ok)  # lost column
  expect_false(tables_equivalence(a, make_tsv(base_df[-2, ]))$ok)  # lost row
  expect_false(tables_equivalence(a, "/nope/missing.txt")$ok)
})

test_that("drop_markers excludes only the named rows and keeps the structural checks", {
  a <- make_tsv(base_df)
  drift <- base_df; drift$BETA[2] <- "0.6"          # rs2 moves by ~17%
  # positive control FIRST: without the exclusion this must fail, or the test
  # below proves nothing
  expect_false(tables_equivalence(a, make_tsv(drift))$ok)
  # dropping the row that moved passes; dropping a DIFFERENT row still fails
  expect_true(tables_equivalence(a, make_tsv(drift), drop_markers = "rs2")$ok)
  expect_false(tables_equivalence(a, make_tsv(drift), drop_markers = "rs3")$ok)
  # null control: dropping nothing is the old behaviour
  expect_true(tables_equivalence(a, make_tsv(base_df), drop_markers = character())$ok)
  # the structural checks survive the exclusion: a table whose MarkerID column
  # differs, or that lacks a named marker, fails rather than silently dropping
  relabelled <- base_df; relabelled$MarkerID[3] <- "rs9"
  expect_false(tables_equivalence(a, make_tsv(relabelled), drop_markers = "rs2")$ok)
  expect_false(tables_equivalence(a, make_tsv(base_df), drop_markers = "rs404")$ok)
  # and a row that differs OUTSIDE the dropped set is still caught even when the
  # dropped row also differs
  both <- drift; both$p.value[1] <- "1.000000E-01"
  expect_false(tables_equivalence(a, make_tsv(both), drop_markers = "rs2")$ok)
})

test_that("comparators work on a real golden file", {
  skip_if_no_goldens()
  g <- file.path(saige_goldens, "step2", "plink_marker_binary.txt")
  cp <- tempfile(fileext = ".txt"); file.copy(g, cp)
  expect_true(files_byte_identical(g, cp))
  expect_true(tables_equivalence(g, cp)$ok)

  # perturb one numeric field in the copy -> both comparators must fail
  lines <- readLines(cp)
  fields <- strsplit(lines[2], "\t")[[1]]
  pcol <- which(strsplit(lines[1], "\t")[[1]] == "p.value")
  fields[pcol] <- format(as.numeric(fields[pcol]) * 1.001, scientific = TRUE)
  lines[2] <- paste(fields, collapse = "\t")
  writeLines(lines, cp)
  expect_false(files_byte_identical(g, cp))
  expect_false(tables_equivalence(g, cp, tol = 1e-5)$ok)
})

test_that("the expect_* wrappers surface pure-layer failures", {
  a <- make_tsv(base_df)
  drift <- base_df; drift$BETA[1] <- "-1.09790"
  b <- make_tsv(drift)
  expect_failure(expect_files_byte_identical(a, b))
  expect_failure(expect_tables_equivalent(a, b, tol = 1e-5))
  expect_success(expect_tables_equivalent(a, make_tsv(base_df)))
})
