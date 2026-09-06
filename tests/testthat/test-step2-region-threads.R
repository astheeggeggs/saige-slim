# --nThreads on a group run threads the C++ region marker loop (2026-09-06).
#
# One marker per worker inside each block of markers_per_chunk, the
# accumulation serial in marker order, so the output is byte-identical to a
# one-thread run whatever the thread count -- one code path over a disjoint
# partition, no cross-marker arithmetic on the workers (the same argument as
# test-step2-threads-inprocess.R). Each run must also REPORT its thread count:
# a knob that never arrived would otherwise pass as a green threaded run
# (LEDGER #30). Five goldens cover the arms that differ on the worker: SKAT-O
# with fastTest escalation (P2), conditional (P3), BURDEN (P5), quantitative
# BURDEN and SKAT-O (P36, P37) -- and the ER arm (MAC <= max_MAC_for_ER on a
# binary trait) is main-thread only and is deferred and re-run there, which P2
# exercises on the bundled fixture.

RT_G <- c("--bedFile=./input/genotype_100markers.bed",
          "--bimFile=./input/genotype_100markers.bim",
          "--famFile=./input/genotype_100markers.fam",
          "--AlleleOrder=alt-first")
RT_SG <- c(
  paste0("--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
         "randomMarkersUsed.sparseGRM.mtx"),
  paste0("--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_",
         "1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
)
RT_BIN <- c("--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
            "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt")
# P36/P37 run WITHOUT the sparse GRM: the quantitative model's variance-ratio
# file has only `null` rows, and --sparseGRMFile needs `sparse` ones.
RT_QUANT <- c(
  paste0("--GMMATmodelFile=", file.path(saige_goldens, "step1", "quant_sparseGRM.rda")),
  paste0("--varianceRatioFile=",
         file.path(saige_goldens, "step1", "quant_sparseGRM.varianceRatio.txt"))
)
RT_GROUP <- c("--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
              "--groupFile=./input/group_new_snpid.txt",
              "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
              "--maxMAF_in_groupTest=0.0001,0.001,0.01")

# Runs the CLI and returns the output path with the log lines attached; the
# expectations stay in the test blocks (lintr cannot see helper-saige.R's
# functions from inside a function body).
cli <- run_saige_cli
cli_ok <- expect_cli_ok
RT_GOLD <- file.path(saige_goldens, "step2")

rt_run <- function(model, extra, nThreads, tag) {
  out <- file.path(tempfile(paste0("step2_rt_", tag, "_", nThreads, "_")), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  log <- tempfile(fileext = ".log")
  st <- cli("step2_SPAtests.R",
            c(RT_G, model, RT_GROUP, extra,
              paste0("--nThreads=", nThreads),
              paste0("--SAIGEOutputFile=", out)),
            log = log)
  structure(out, status = st, log = readLines(log, warn = FALSE), nThreads = nThreads)
}

# Every threaded run must say so: a knob that never arrived would otherwise pass
# as a green threaded run (LEDGER #30).
expect_rt_ok <- function(out) {
  cli_ok(attr(out, "status"))
  wanted <- paste0("Region marker loop will use ", attr(out, "nThreads"),
                   " in-process threads")
  testthat::expect_true(any(grepl(wanted, attr(out, "log"), fixed = TRUE)))
}

test_that("SKAT-O with fastTest escalation is byte-identical on 4 threads (P2)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- rt_run(RT_BIN, c(RT_SG, "--is_fastTest=TRUE"), 4, "p2")
  expect_rt_ok(out)
  expect_output_matches_golden(file.path(RT_GOLD, "plink_group_fastTest.txt"), out)
})

test_that("the conditional region path is byte-identical on 4 threads (P3)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- rt_run(RT_BIN, c(RT_SG, "--condition=rs30,rs79", "--is_fastTest=TRUE"), 4, "p3")
  expect_rt_ok(out)
  expect_output_matches_golden(file.path(RT_GOLD, "plink_group_cond_binary.txt"), out)
})

test_that("BURDEN is byte-identical on 8 threads (P5)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- rt_run(RT_BIN, c(RT_SG, "--r.corr=1", "--is_fastTest=FALSE"), 8, "p5")
  expect_rt_ok(out)
  expect_output_matches_golden(file.path(RT_GOLD, "plink_group_burden.txt"), out)
})

test_that("the quantitative BURDEN and SKAT-O arms are byte-identical on 4 threads (P36, P37)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- rt_run(RT_QUANT, c("--r.corr=1", "--is_fastTest=FALSE"), 4, "p36")
  expect_rt_ok(out)
  expect_output_matches_golden(file.path(RT_GOLD, "plink_group_burden_quant.txt"), out)
  out <- rt_run(RT_QUANT, "--is_fastTest=FALSE", 4, "p37")
  expect_rt_ok(out)
  expect_output_matches_golden(file.path(RT_GOLD, "plink_group_nofast_quant.txt"), out)
})
