# LOCO end to end, which no golden reached before 2026-09-08.
#
# `--LOCO=TRUE` appeared in 0 of the 51 golden runs and `--LOCO=FALSE` in 46.
# LOCO is three mechanisms (CLAUDE.md): step-1 global-state mutation per
# chromosome, the save-time file layout, and step-2 `ReadModel` overwriting the
# model's fitted values from `LOCOResult[[chrom]]`. `test-step1-lowmem-loco.R`
# covers the second and `test-unit-readmodel.R` the third at unit level; this
# file covers the third END TO END, against golden `binary_loco_2chr`.
#
# Only a dense-GRM fit can carry a LOCOResult at all:
# `useSparseGRMtoFitNULL=TRUE` forces `LOCO <- FALSE`
# (R/SAIGE_fitGLMM_fast.R:836).

LOCO_G2 <- function() {
  c("--bedFile=./input/genotype_100markers_2chr.bed",
    "--bimFile=./input/genotype_100markers_2chr.bim",
    "--famFile=./input/genotype_100markers_2chr.fam",
    "--AlleleOrder=alt-first")
}
loco_model <- function() {
  c(paste0("--GMMATmodelFile=",
           file.path(saige_goldens, "step1", "binary_loco_2chr.rda")),
    paste0("--varianceRatioFile=",
           file.path(saige_goldens, "step1", "binary_loco_2chr.varianceRatio.txt")))
}
LOCO_COMMON <- c("--minMAF=0", "--minMAC=1", "--is_output_moreDetails=TRUE")

loco_cache <- new.env(parent = emptyenv())
loco_run <- function(key, extra) {
  if (!is.null(loco_cache[[key]])) return(loco_cache[[key]])
  out <- file.path(tempfile("loco_"), paste0(key, ".txt"))
  dir.create(dirname(out), recursive = TRUE)
  status <- run_saige_cli("step2_SPAtests.R",
                          c(LOCO_G2(), loco_model(), LOCO_COMMON, extra,
                            paste0("--SAIGEOutputFile=", out)))
  expect_cli_ok(status)
  loco_cache[[key]] <- read_saige_table(out)
  loco_cache[[key]]
}

test_that("the golden LOCO model carries a per-chromosome refit, and only where markers are", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  e <- new.env()
  load(file.path(saige_goldens, "step1", "binary_loco_2chr.rda"), envir = e)
  m <- e$modglmm
  expect_true(m$LOCO)
  expect_length(m$LOCOResult, 22L)
  flag <- vapply(m$LOCOResult, function(x) isTRUE(x$isLOCO), logical(1))
  # The fixture has markers on chromosomes 1 and 2 only; the other 20 entries
  # are list(isLOCO = FALSE) stubs, which is what keeps the .rda at 677 KB.
  expect_equal(which(flag), 1:2)
  a <- m$LOCOResult[[1]]
  b <- m$LOCOResult[[2]]
  # The two refits must actually differ, or the whole mechanism is a no-op and
  # every assertion below would pass against a model that ignores --chrom.
  expect_false(isTRUE(all.equal(a$coefficients, b$coefficients)))
  expect_gt(max(abs(a$residuals - b$residuals)), 1e-3)
})

test_that("step 2 applies LOCOResult[[chrom]]: LOCO and non-LOCO differ on every marker", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # Same model, same chromosome, only --LOCO changed. If ReadModel did not swap
  # the fitted values in, these would be identical.
  for (chrom in c("1", "2")) {
    on_  <- loco_run(paste0("on", chrom),  c(paste0("--chrom=", chrom), "--LOCO=TRUE"))
    off_ <- loco_run(paste0("off", chrom), c(paste0("--chrom=", chrom), "--LOCO=FALSE"))
    expect_equal(on_$MarkerID, off_$MarkerID)
    expect_true(all(as.numeric(on_$p.value) != as.numeric(off_$p.value)),
                label = paste0("every p.value moves on chrom ", chrom))
  }
})

test_that("leaving out the marker-dense chromosome moves p-values far more", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  relmove <- function(chrom) {
    on_  <- loco_run(paste0("on", chrom),  c(paste0("--chrom=", chrom), "--LOCO=TRUE"))
    off_ <- loco_run(paste0("off", chrom), c(paste0("--chrom=", chrom), "--LOCO=FALSE"))
    median(abs(as.numeric(on_$p.value) - as.numeric(off_$p.value)) /
             pmax(abs(as.numeric(off_$p.value)), 1e-300))
  }
  # The GRM fixture is 9,988 markers on chr1 against 12 on chr2, so leaving out
  # chr1 changes the model enormously and leaving out chr2 barely at all. This
  # is a claim about the DIRECTION of the asymmetry, predicted from the marker
  # counts before it was measured, not about the exact ratio.
  expect_gt(relmove("1"), relmove("2"))
})

test_that("--chrom is required under LOCO, and the message names the wrong test type", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- file.path(tempfile("loco_nochrom_"), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  log <- tempfile(fileext = ".log")
  status <- run_saige_cli("step2_SPAtests.R",
                          c(LOCO_G2(), loco_model(), LOCO_COMMON, "--LOCO=TRUE",
                            paste0("--SAIGEOutputFile=", out)), log = log)
  expect_false(status == 0)
  txt <- paste(readLines(log, warn = FALSE), collapse = "\n")
  # Characterization: this is a SINGLE-VARIANT run and the message says
  # "gene- or region-based tests". Pinned as-is; see LEDGER.
  expect_match(txt, "chrom needs to be specified", fixed = TRUE)
  expect_match(txt, "gene- or region-based tests", fixed = TRUE)
})

test_that("a non-LOCO model under --LOCO=TRUE is refused, and that message IS accurate", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # The sibling guard three lines above the one in the previous test, on the
  # same path. It is pinned too so that #110 records one bad string rather than
  # leaving the reader to wonder whether the whole block is wrong.
  out <- file.path(tempfile("loco_nonloco_"), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  log <- tempfile(fileext = ".log")
  status <- run_saige_cli("step2_SPAtests.R",
                          c(LOCO_G2(), LOCO_COMMON, "--chrom=1", "--LOCO=TRUE",
                            paste0("--GMMATmodelFile=",
                                   file.path(saige_goldens, "step1", "binary_denseGRM.rda")),
                            paste0("--varianceRatioFile=",
                                   file.path(saige_goldens, "step1",
                                             "binary_denseGRM.varianceRatio.txt")),
                            paste0("--SAIGEOutputFile=", out)), log = log)
  expect_false(status == 0)
  expect_match(paste(readLines(log, warn = FALSE), collapse = "\n"),
               "does not contain LOCO results", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# LEDGER #111 and #112: the LOCO defaults, and LOCO against a sparse-GRM model.
# Characterization -- these assert what the code does TODAY, including the
# contradiction and the unfollowable advice. Invert them when either is fixed.

test_that("step 2's --LOCO help contradicts its own default (LEDGER #111)", {
  skip_unless_level("fast")
  cli <- file.path(saige_extdata, "step2_SPAtests.R")
  skip_if_not(file.exists(cli))
  src <- paste(readLines(cli, warn = FALSE), collapse = "\n")
  # The declaration and the help text sit four lines apart in one make_option.
  block <- regmatches(src, regexpr('"--LOCO".*?\\)\\s*,', src))
  expect_match(block, "default = TRUE", fixed = TRUE)
  expect_match(block, "[default=FALSE]", fixed = TRUE)
})

test_that("LOCO against a sparse-GRM model is refused, with advice it cannot follow (#112)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # binary_sparseGRM_a was fitted with --useSparseGRMtoFitNULL=TRUE, which sets
  # LOCO <- FALSE silently (R/SAIGE_fitGLMM_fast.R:836) even though that run
  # requested --LOCO=TRUE. Here --LOCO is simply left at its CLI default, which
  # IS TRUE (#111), and --chrom is supplied so this is not #110's path. So this
  # is the default invocation of the recommended workflow.
  out <- file.path(tempfile("loco_sparse_"), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  log <- tempfile(fileext = ".log")
  args <- c("--bedFile=./input/genotype_100markers.bed",
            "--bimFile=./input/genotype_100markers.bim",
            "--famFile=./input/genotype_100markers.fam",
            "--AlleleOrder=alt-first",
            paste0("--GMMATmodelFile=",
                   file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
            paste0("--varianceRatioFile=",
                   file.path(saige_goldens, "step1",
                             "binary_sparseGRM_a.varianceRatio.txt")),
            "--chrom=1", "--minMAF=0", "--minMAC=1",
            "--is_output_moreDetails=TRUE",
            paste0("--SAIGEOutputFile=", out))
  status <- run_saige_cli("step2_SPAtests.R", args, log = log)
  expect_false(status == 0)
  txt <- paste(readLines(log, warn = FALSE), collapse = "\n")
  expect_match(txt, "does not contain LOCO results", fixed = TRUE)
  # The half of the remedy that cannot be carried out on a sparse-GRM fit: the
  # user DID run step 1 with LOCO and it was overridden without a word.
  expect_match(txt, "please run Step 1 using LOCO", fixed = TRUE)
  # ... and the half that can.
  expect_match(txt, "please set LOCO=FALSE in this step", fixed = TRUE)
})
