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

# Golden P6's command (group test, SKAT-O, is_fastTest=FALSE) plus --nThreads.
skato_run <- function(nThreads, extra = character()) {
  out <- file.path(tempfile(paste0("step2_skato_", nThreads, "_")), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- cli2("step2_SPAtests.R", c(
    SKATO_G, SKATO_MODEL, SKATO_SG, extra,
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=FALSE",
    paste0("--nThreads=", nThreads),
    paste0("--SAIGEOutputFile=", out)
  ))
  list(status = st, path = out,
       log = readLines(attr(st, "log"), warn = FALSE))
}

workers_reported <- function(lg) {
  sum(grepl("SKAT-O post-processing will use", lg, fixed = TRUE))
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
  # And it says plainly that the marker loop is not what got threaded.
  expect_true(any(grepl("does not thread the group test's marker loop",
                        r$log, fixed = TRUE)))
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
  expect_output_matches_golden(P3g, r$path)
})
