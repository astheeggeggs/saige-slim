# LEDGER #52: the psam sample map keyed on column 0, which is FID.
#
# `PgenShared::readPsamFile` did `sampleInPgen[line_elements[0]] = idx;  // IID`.
# Column 0 of a "#FID IID SEX" psam is FID, so step 2 failed outright on any
# psam with FID != IID, and two samples sharing an FID collapsed into one map
# entry -- which also undercounts N0, since that is `sampleInPgen.size()`.
#
# EVERY committed PGEN fixture has FID == IID, which is the only reason all 13
# PGEN step-2 goldens ever passed. The gate was the missing fixture, not the
# fix: with FID == IID the bug is invisible no matter how many goldens exist.
#
# So the fixtures are DERIVED here rather than committed, which also keeps the
# FID/IID relationship under test visible in the source instead of buried in a
# 10,000-line file. The pgen and pvar are the committed ones; only the psam is
# rewritten, and rewriting it must not change a single output byte, because the
# genotypes and the sample ORDER are untouched.
#
# Tier: standard.

PSAM_PGV <- c("--pgenFile=./input/genotype_100markers.pgen",
              "--pvarFile=./input/genotype_100markers.pvar")

# Aliases so the helpers below reference names this file defines. `run_saige_cli`
# and the two path constants live in helper-saige.R, which lintr does not see.
cli <- run_saige_cli
goldensDir <- saige_goldens
extdataDir <- saige_extdata

P13 <- file.path(goldensDir, "step2", "pgen_marker_binary.txt")
BIN <- c(
  paste0("--GMMATmodelFile=",
         file.path(goldensDir, "step1", "binary_sparseGRM_a.rda")),
  paste0("--varianceRatioFile=",
         file.path(goldensDir, "step1", "binary_sparseGRM_a.varianceRatio.txt"))
)

# The committed psam, as a header vector plus a character matrix of data rows.
psam_cols <- function() {
  ln <- readLines(file.path(extdataDir, "input", "genotype_100markers.psam"))
  hdr <- strsplit(sub("^#", "", ln[1]), "\t")[[1]]
  body <- do.call(rbind, strsplit(ln[-1], "\t"))
  colnames(body) <- hdr
  list(hdr = hdr, body = body)
}

write_psam <- function(hdr, body) {
  out <- tempfile(fileext = ".psam")
  writeLines(c(paste0("#", paste(hdr, collapse = "\t")),
               apply(body, 1, paste, collapse = "\t")), out)
  out
}

# One marker scan (golden P13's command) against a caller-supplied psam.
run_with_psam <- function(psam, tag) {
  out <- file.path(tempfile(paste0("step2_psam_", tag, "_")), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- cli("step2_SPAtests.R", c(
    PSAM_PGV, paste0("--psamFile=", psam),
    "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    BIN, "--is_output_moreDetails=TRUE",
    paste0("--SAIGEOutputFile=", out)
  ))
  list(status = st, path = out)
}

test_that("the committed psam still reproduces golden P13 when rewritten verbatim", {
  # Negative control for the two tests below: the rewrite machinery itself
  # (readLines -> split -> paste -> writeLines) must be output-neutral, or a
  # pass on the FID != IID fixture would say nothing about the FID != IID part.
  skip_unless_level("standard")
  skip_if_no_goldens()

  p <- psam_cols()
  r <- run_with_psam(write_psam(p$hdr, p$body), "verbatim")
  expect_cli_ok(r$status)
  expect_output_matches_golden(P13, r$path)
})

test_that("a psam whose FID differs from its IID gives golden P13", {
  skip_unless_level("standard")
  skip_if_no_goldens()

  p <- psam_cols()
  body <- p$body
  # One FID shared by all 10,000 samples: FID != IID for every row AND every
  # row collides on column 0, which is the undercount half of #52. At unpatched
  # HEAD the map holds ONE entry and the run dies on the first sample lookup.
  body[, "FID"] <- "FAM1"
  expect_false(any(body[, "FID"] == body[, "IID"]))

  r <- run_with_psam(write_psam(p$hdr, body), "fidne")
  expect_cli_ok(r$status)
  expect_output_matches_golden(P13, r$path)
})

test_that("a psam with no FID column at all gives golden P13", {
  skip_unless_level("standard")
  skip_if_no_goldens()

  # "#IID SEX" is what plink2 writes when there are no family IDs. Here IID IS
  # column 0, so this passes at unpatched HEAD too -- it is here to pin that
  # finding the column BY NAME did not break the no-FID layout.
  p <- psam_cols()
  keep <- setdiff(p$hdr, "FID")
  r <- run_with_psam(write_psam(keep, p$body[, keep, drop = FALSE]), "nofid")
  expect_cli_ok(r$status)
  expect_output_matches_golden(P13, r$path)
})

test_that("a psam with a duplicate IID is refused rather than silently undercounted", {
  skip_unless_level("standard")

  p <- psam_cols()
  body <- p$body
  body[2, "IID"] <- body[1, "IID"]
  r <- run_with_psam(write_psam(p$hdr, body), "dupiid")
  expect_false(identical(as.integer(r$status), 0L))
  lg <- readLines(attr(r$status, "log"), warn = FALSE)
  expect_true(any(grepl("duplicate IID", lg, fixed = TRUE)))
})

test_that("a psam with no column header is refused", {
  skip_unless_level("standard")

  p <- psam_cols()
  out <- tempfile(fileext = ".psam")
  writeLines(apply(p$body, 1, paste, collapse = "\t"), out)
  r <- run_with_psam(out, "nohdr")
  expect_false(identical(as.integer(r$status), 0L))
  lg <- readLines(attr(r$status, "log"), warn = FALSE)
  expect_true(any(grepl("column header", lg, fixed = TRUE)))
})
