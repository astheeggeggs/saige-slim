# LEDGER #155 (FIXED 2026-09-15): the SKAT-O pool retries its port.
#
# parallel::makeForkCluster listens on ONE port fixed at package load from the
# RNG and the clock (11000-11999). Two SAIGE processes on one box can draw the
# same one, and the second died in startSkatoCluster with `creation of server
# socket failed: port NNNNN cannot be opened` before any test -- 1 in 200 of
# the #154 stress drill's concurrent runs, eight at a time, which is the
# pre-push gate's own shape. Now a failed bind is retried on a port derived
# from the pid and the attempt number.
#
# Characterized by holding a port open in a helper process and forcing the
# run onto it with R_PARALLEL_PORT (parallel reads it at load). Before the
# fix: the message above at exit 1. Now: the retry line, exit 0, results
# byte-identical to golden P3. Tier: standard (one 4-thread region run).

l155_cli <- function(...) get("run_saige_cli")(...)

test_that("a taken cluster port is retried, not fatal (LEDGER #155)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  port <- 11000L + as.integer(Sys.getpid() %% 1000)
  # a helper R process that holds the port for up to two minutes
  rscript <- file.path(R.home("bin"), "Rscript")
  holder <- tempfile("l155_holder_", fileext = ".log")
  pid_file <- tempfile("l155_pid_")
  system2(rscript, c("-e", shQuote(sprintf(
    "writeLines(as.character(Sys.getpid()), '%s'); s <- serverSocket(%dL); Sys.sleep(120)",
    pid_file, port))), stdout = holder, stderr = holder, wait = FALSE)
  for (i in 1:50) { if (file.exists(pid_file)) break; Sys.sleep(0.1) }
  Sys.sleep(0.5)
  on.exit(try(tools::pskill(as.integer(readLines(pid_file, warn = FALSE))), silent = TRUE),
          add = TRUE)
  out <- file.path(tempfile("step2_l155_"), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- l155_cli("step2_SPAtests.R", c(
    "--bedFile=./input/genotype_100markers.bed", "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
    "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
    "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
    "--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
    paste0("--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01",
    "--condition=rs30,rs79", "--is_fastTest=TRUE", "--nThreads=4",
    paste0("--SAIGEOutputFile=", out)
  ), env = c(R_PARALLEL_PORT = as.character(port)))
  lg <- readLines(attr(st, "log"), warn = FALSE)
  # the control that the port really was taken: the first bind failed
  expect_true(any(grepl(paste0("port ", port, " cannot be opened"), lg, fixed = TRUE)))
  expect_true(any(grepl("retrying on port", lg, fixed = TRUE)))
  expect_identical(as.integer(st), 0L)
  expect_true(any(grepl("Analysis done!", lg, fixed = TRUE)))
  expect_output_matches_golden(file.path(saige_goldens, "step2", "plink_group_cond_binary.txt"), out)
})
