# Parallel SKAT-O post-processing (--nThreads on a group run).
#
# `Rprof` at N = 50,000 puts 45% of a group run in `SKAT_davies` under
# `integrate()` -- R adaptive quadrature over the SKAT package's compiled Davies
# routine. Each (annotation x max_MAF x weight set) cell's SKAT-O p-value is a
# pure function of that cell's (Score, Phi, r.corr): no cross-cell reduction, no
# shared accumulator, no RNG. So this is ONE code path re-run over a partition,
# which is the narrow case where byte-identity is the right gate rather than a
# category error (cf. CLAUDE.md's rule against comparing two different paths).
#
# The tests below assert byte-identity AND that the parallel path was actually
# taken. Without the second, a run whose option never arrived -- or which fell
# through `length(st$jobs) < 2L` for every region -- would pass as a green
# parallel run while proving nothing. That is the LEDGER #30 shape and it is
# what `test-step2-threads-inprocess.R` guards against on the marker side.

SKATO_G <- c("--bedFile=./input/genotype_100markers.bed",
             "--bimFile=./input/genotype_100markers.bim",
             "--famFile=./input/genotype_100markers.fam",
             "--AlleleOrder=alt-first")
SKATO_MODEL <- c(
  "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
  "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt"
)
SKATO_SG <- c(
  paste0("--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
         "randomMarkersUsed.sparseGRM.mtx"),
  paste0("--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_",
         "1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
)
cli2 <- run_saige_cli
P6g <- file.path(saige_goldens, "step2", "plink_group_nofast.txt")
P3g <- file.path(saige_goldens, "step2", "plink_group_cond_binary.txt")
P35g <- file.path(saige_goldens, "step2", "plink_group_fastTest_escalated.txt")

# Golden P6's command (group test, SKAT-O, is_fastTest=FALSE) plus --nThreads.
# `extra` OVERRIDES a default of the same name: optparse keeps the FIRST value
# of a repeated option, so appending "--is_fastTest=TRUE" after the default
# "--is_fastTest=FALSE" used to leave the run on the non-fast path while the
# test read as a fastTest one (P2/P35's bytes equal P6's only where LEDGER #51
# says so, which is why it passed).
skato_run <- function(nThreads, extra = character(), env = character()) {
  out <- file.path(tempfile(paste0("step2_skato_", nThreads, "_")), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  defaults <- c(
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=FALSE"
  )
  opt_name <- function(a) sub("=.*$", "", a)
  defaults <- defaults[!opt_name(defaults) %in% opt_name(extra)]
  st <- withr::with_envvar(env, cli2("step2_SPAtests.R", c(
    SKATO_G, SKATO_MODEL, SKATO_SG, extra, defaults,
    paste0("--nThreads=", nThreads),
    paste0("--SAIGEOutputFile=", out)
  )))
  list(status = st, path = out,
       log = readLines(attr(st, "log"), warn = FALSE))
}

workers_reported <- function(lg) {
  sum(grepl("SKAT-O post-processing will use", lg, fixed = TRUE))
}
# The batch lines (step 5, 2026-09-06): "SKAT-O batch: N cell(s) from M
# region(s) to the pool", one per batch; regions per batch as integers.
regions_per_batch <- function(lg) {
  b <- grep("^SKAT-O batch: .* from [0-9]+ region", lg, value = TRUE)
  as.integer(sub("^SKAT-O batch: .* from ([0-9]+) region.*$", "\\1", b))
}

test_that("nThreads=1 leaves the group path untouched", {
  skip_unless_level("standard")
  skip_if_no_goldens()

  r <- skato_run(1)
  expect_cli_ok(r$status)
  # The pool is never built, so the serial code path is literally unchanged.
  expect_identical(workers_reported(r$log), 0L)
  expect_output_matches_golden(P6g, r$path)
})

test_that("nThreads>1 evaluates SKAT-O in parallel and is byte-identical", {
  skip_unless_level("standard")
  skip_if_no_goldens()

  r <- skato_run(4)
  expect_cli_ok(r$status)
  # It really ran in parallel: the pool is announced once, on creation, and
  # only a region that took the parallel branch creates it.
  expect_identical(workers_reported(r$log), 1L)
  expect_true(any(grepl("use 4 workers", r$log, fixed = TRUE)))
  # And the region marker loop threads too, and says so (2026-09-06).
  expect_true(any(grepl("Region marker loop will use 4 in-process threads",
                        r$log, fixed = TRUE)))
  # Step 5: the cells of the whole chunk went to the pool in ONE batch. The
  # bundled group file has two regions in one chunk; a batch of one region
  # would be the 2026-09-04 per-region pool wearing the new label.
  expect_identical(regions_per_batch(r$log), 2L)
  expect_output_matches_golden(P6g, r$path)
})

test_that("a zero batch budget falls back to one region per batch, byte-identical", {
  skip_unless_level("standard")
  skip_if_no_goldens()

  # The other way round: with the memory budget at zero every region is its
  # own batch, so the two regions produce two batch lines -- and the same
  # bytes, since the batch boundary cannot reach a cell's answer.
  r <- skato_run(4, env = c(SAIGE_SKATO_BATCH_MB = "0"))
  expect_cli_ok(r$status)
  expect_identical(workers_reported(r$log), 1L)
  expect_identical(regions_per_batch(r$log), c(1L, 1L))
  expect_output_matches_golden(P6g, r$path)
})

test_that("the conditional group path is byte-identical in parallel too", {
  skip_unless_level("standard")
  skip_if_no_goldens()

  # --condition takes the SECOND skat_fun call site, inside the isCondition
  # branch, which the unconditional runs above never reach. Golden P3.
  extra <- c("--condition=rs30,rs79", "--is_fastTest=TRUE")
  r <- skato_run(4, extra)
  expect_cli_ok(r$status)
  expect_identical(workers_reported(r$log), 1L)
  expect_identical(regions_per_batch(r$log), 2L)
  expect_output_matches_golden(P3g, r$path)
})

test_that("fastTest escalation inside a batch keeps the single rows in region order (P35)", {
  skip_unless_level("standard")
  skip_if_no_goldens()

  # Golden P35: SKAT-O with fastTest and the cutoff at 1, so EVERY region's
  # Cauchy p crosses it and the accurate pass re-runs each (at the default
  # cutoff no bundled gene escalates, LEDGER #51). Batched, the decision comes
  # after the batch's pool call, so every region's fast-pass single-variant
  # rows are stashed and the escalated models' cells go round the pool once
  # more. P35's singleAssoc file pins the OTHER writer's rows (LEDGER #50),
  # the region file the p-values.
  r <- skato_run(4, c("--is_fastTest=TRUE", "--pval_cutoff_for_fastTest=1"))
  expect_cli_ok(r$status)
  expect_identical(workers_reported(r$log), 1L)
  expect_identical(regions_per_batch(r$log), 2L)
  expect_identical(sum(grepl("Non-fast test is performed", r$log, fixed = TRUE)), 2L)
  expect_true(any(grepl("escalated cell(s) to the pool", r$log, fixed = TRUE)))
  expect_output_matches_golden(P35g, r$path)
  expect_output_matches_golden(paste0(P35g, ".singleAssoc.txt"),
                               paste0(r$path, ".singleAssoc.txt"))
})
