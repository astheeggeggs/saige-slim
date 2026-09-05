# LEDGER #88 and #89 FIXED: the variant-file parsers refuse malformed input by
# file, line and column, and stop refusing input they never used.
#
# #88: readPvarFile had none of readPsamFile's hardening. A blank line (one
# newline appended to the committed pvar) made it index a split vector out of
# range and then die in std::stoi with "stoi: no conversion" -- naming nothing.
# It now skips blank lines, strips \r, refuses a short row, and refuses a data
# row before the #CHROM header.
# #89: three bare std::sto* calls. bim column 3 (centimorgans) was parsed into
# a member nothing read, so '.' there refused a usable dataset; it is not
# parsed any more. bim column 4 and pvar POS are parsed with a message naming
# file, line and column. Fixtures are DERIVED from the committed files.
# Tier: standard.

vp_m <- function() file.path(saige_goldens, "step1", "binary_sparseGRM_a")
vp_common <- function(out) {
  m <- vp_m()
  c(paste0("--SAIGEOutputFile=", out), "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    paste0("--GMMATmodelFile=", m, ".rda"), paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
    "--is_output_moreDetails=TRUE")
}
vp_pgen <- function(pvar) {
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", c("--pgenFile=./input/genotype_100markers.pgen",
    paste0("--pvarFile=", pvar), "--psamFile=./input/genotype_100markers.psam", vp_common(out)))
  list(status = as.integer(st), log = readLines(attr(st, "log"), warn = FALSE), out = out)
}
vp_plink <- function(bim) {
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", c("--bedFile=./input/genotype_100markers.bed",
    paste0("--bimFile=", bim), "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first", vp_common(out)))
  list(status = as.integer(st), log = readLines(attr(st, "log"), warn = FALSE), out = out)
}
derive <- function(src, ext, f) {
  lines <- readLines(src, warn = FALSE)
  out <- tempfile(fileext = ext); writeLines(f(lines), out); out
}
pvar_src <- function() file.path(saige_extdata, "input", "genotype_100markers.pvar")
bim_src <- function() file.path(saige_extdata, "input", "genotype_100markers.bim")

test_that("#88: a blank line in the pvar is skipped and the run matches the PGEN golden", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- vp_pgen(derive(pvar_src(), ".pvar", function(l) c(l, "")))
  expect_identical(r$status, 0L)
  expect_output_matches_golden(file.path(saige_goldens, "step2", "pgen_marker_binary.txt"), r$out)
})

test_that("#88: a short pvar row is refused by file and line", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  short <- function(l) { i <- grep("^[^#]", l)[5]; l[i] <- "1\t5"; l }   # data row 5 cut to two fields
  r <- vp_pgen(derive(pvar_src(), ".pvar", short))
  expect_false(identical(r$status, 0L))
  expect_true(any(grepl("pvar file line [0-9]+ has 2 fields but the header needs at least", r$log)))
  expect_true(any(grepl("genotype_100markers.pvar|\\.pvar\\)", r$log)))
})

test_that("#89: a non-numeric pvar POS is refused naming file, line and column", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  bad <- function(l) { i <- grep("^[^#]", l)[7]; f <- strsplit(l[i], "\t")[[1]]; f[2] <- "12abc"; l[i] <- paste(f, collapse = "\t"); l }
  r <- vp_pgen(derive(pvar_src(), ".pvar", bad))
  expect_false(identical(r$status, 0L))
  expect_true(any(grepl("Cannot parse the POS field '12abc' on line [0-9]+ of", r$log)))
  expect_false(any(grepl("stoi", r$log, fixed = TRUE)))
})

test_that("#89: '.' in bim column 3 (centimorgans, never read) no longer refuses the dataset", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  dot <- function(l) vapply(l, function(x) { f <- strsplit(x, "\t")[[1]]; f[3] <- "."; paste(f, collapse = "\t") }, "", USE.NAMES = FALSE)
  r <- vp_plink(derive(bim_src(), ".bim", dot))
  expect_identical(r$status, 0L)
  expect_output_matches_golden(file.path(saige_goldens, "step2", "plink_marker_binary.txt"), r$out)
})

test_that("#89: a non-numeric bim position, and a short bim row, are refused by file and line", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  bad <- function(l) { f <- strsplit(l[9], "\t")[[1]]; f[4] <- "9x"; l[9] <- paste(f, collapse = "\t"); l }
  r <- vp_plink(derive(bim_src(), ".bim", bad))
  expect_false(identical(r$status, 0L))
  expect_true(any(grepl("Cannot parse the position field '9x' on line 9 of .*column 4 of a .bim", r$log)))
  short <- function(l) { l[4] <- "1\trs4"; l }
  r2 <- vp_plink(derive(bim_src(), ".bim", short))
  expect_false(identical(r2$status, 0L))
  expect_true(any(grepl("bim file line 4 has 2 fields; a .bim row needs 6", r2$log, fixed = TRUE)))
})
