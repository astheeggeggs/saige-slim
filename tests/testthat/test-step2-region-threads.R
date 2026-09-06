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

# Step 3b (2026-09-06) keeps the per-cell test records of the ultra-rare pass
# and the BURDEN emitter in storage that persists across regions, so a region's
# output must not depend on the regions before it. The first persistent build
# re-emitted an earlier region's result for a cell the current region did not
# test (its `tested` flag was never reset); the bundled BURDEN golden did not
# see it because both its regions test the same cells, and the sim50k byte gate
# did. Region B here is the bundled GENE2's lof markers plus two of its
# missense markers, rs55 and rs62, each with ONE minor allele among the 1,000
# analysed samples: under --annotation_in_groupTest=lof,missense B's missense
# cell is present (so R keeps the annotation and the cell count matches A's)
# but untested (MAC 2 < --minGroupMAC_in_BurdenTest 5), right after a region
# A (GENE1) whose missense cell was tested (MAC 5). A missense annotation with
# NO markers would not do: R drops absent annotations per region, the record
# count changes, and the records are rebuilt.
rt_group_lines <- function() {
  g <- readLines("./input/group_new_snpid.txt")
  tok <- strsplit(g, " ", fixed = TRUE)
  gene <- vapply(tok, `[`, "", 1)
  kind <- vapply(tok, `[`, "", 2)
  var <- tok[[which(gene == "GENE2" & kind == "var")]]
  anno <- tok[[which(gene == "GENE2" & kind == "anno")]]
  keep <- c(TRUE, TRUE, anno[-(1:2)] == "lof" | var[-(1:2)] %in% c("rs55", "rs62"))
  list(a = g[gene == "GENE1"],
       b = c(paste(var[keep], collapse = " "), paste(anno[keep], collapse = " ")))
}

rt_run_burden_gf <- function(groupFile, nThreads) {
  out <- file.path(tempfile(paste0("step2_rt_gf_", nThreads, "_")), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  log <- tempfile(fileext = ".log")
  st <- cli("step2_SPAtests.R",
            c(RT_G, RT_BIN, RT_SG,
              "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
              paste0("--groupFile=", groupFile),
              "--annotation_in_groupTest=lof,missense",
              # 0.001 first: region A's first missense cell must be TESTED for
              # a stale flag to exist (no bundled marker has MAF <= 0.0001)
              "--maxMAF_in_groupTest=0.001,0.01",
              "--r.corr=1", "--is_fastTest=FALSE",
              paste0("--nThreads=", nThreads),
              paste0("--SAIGEOutputFile=", out)),
            log = log)
  structure(out, status = st, log = readLines(log, warn = FALSE), nThreads = nThreads)
}

test_that("a region's BURDEN output does not depend on the regions before it (1 and 4 threads)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  gl <- withr::with_dir(saige_extdata, rt_group_lines())
  gf_ab <- tempfile("group_ab_", fileext = ".txt")
  gf_b <- tempfile("group_b_", fileext = ".txt")
  writeLines(c(gl$a, gl$b), gf_ab)
  writeLines(gl$b, gf_b)
  region_b <- function(out) {
    rows <- readLines(out)
    rows[startsWith(rows, "GENE2\t")]
  }
  for (nt in c(1L, 4L)) {
    out_ab <- rt_run_burden_gf(gf_ab, nt)
    out_b <- rt_run_burden_gf(gf_b, nt)
    cli_ok(attr(out_ab, "status"))
    cli_ok(attr(out_b, "status"))
    b_after_a <- region_b(out_ab)
    b_alone <- region_b(out_b)
    testthat::expect_gt(length(b_alone), 1L)
    # region B's missense cell has MAC 2 < 5, so no missense row may be emitted
    testthat::expect_false(any(grepl("\tmissense\t", b_alone, fixed = TRUE)))
    testthat::expect_false(any(grepl("\tmissense\t", b_after_a, fixed = TRUE)))
    testthat::expect_identical(b_after_a, b_alone)
  }
})

# The ER resampling arm on a COLLAPSED cell (a binary pseudo-marker with
# MAC <= max_MAC_for_ER whose exact test is too large to enumerate) needs R's
# RNG, so a worker must defer it to the main thread. No bundled golden reaches
# it (drilled 2026-09-06: a worker that refused instead of deferring left P2 at
# 4 threads and P5 at 8 green), so this run raises --max_MAC_for_ER to 200,
# the value the ER resampling arm needs on this fixture (P33), and holds
# 4 threads byte-identical to 1 for the SKAT-O (pass B) and BURDEN (pass D)
# cell tests. Same command, same path, disjoint partition: byte-identity is
# the right gate.
test_that("collapsed cells on the ER resampling arm are byte-identical on 4 threads", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  arms <- list(
    skato = c(RT_SG, "--max_MAC_for_ER=200", "--is_fastTest=FALSE"),
    burden = c(RT_SG, "--max_MAC_for_ER=200", "--r.corr=1", "--is_fastTest=FALSE")
  )
  for (kind in names(arms)) {
    out1 <- rt_run(RT_BIN, arms[[kind]], 1, paste0("er_", kind))
    out4 <- rt_run(RT_BIN, arms[[kind]], 4, paste0("er_", kind))
    expect_rt_ok(out1)
    expect_rt_ok(out4)
    testthat::expect_identical(readLines(out4), readLines(out1))
    sa1 <- paste0(out1, ".singleAssoc.txt")
    if (file.exists(sa1)) {
      testthat::expect_identical(readLines(paste0(out4, ".singleAssoc.txt")), readLines(sa1))
    }
  }
})
