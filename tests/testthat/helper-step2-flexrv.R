# Shared scaffolding for the FlexRV step-2 tests (test-step2-flexrv*.R).
# Nothing here runs at source time: fx_run() calls skip_unless_level("standard")
# before it spawns anything, so the fast tier pays only the definitions.
#
# The fixture, extdata/input/group_new_snpid_flexrv.txt, is
# group_new_snpid.txt (GENE1: 30 lof + 20 missense, GENE2: 20 lof + 30
# missense) plus one `score:REVEL` line per gene: LoF markers scored 1 (the
# paper's rule), missense markers carrying the REVEL values of
# group_new_snpid_multiWeights.txt, interior values in (0, 1) so every
# nonlinear score transform is live.
#
# The fixture ARITHMETIC (2026-09-09, measured with plink2 --freq on the model's
# 1,000 samples -- NOT the fam file's 10,000; step 2 subsets to the model):
# 2N = 2,000, so max MAF 0.01 admits MAC <= 20. At that cell GENE1 holds 21
# polymorphic markers with MAC 1..20 and GENE2 18 with MAC 1..13; at the
# default collapse cutoff (10) the real markers are rs4/rs17/rs31/rs33 (GENE1)
# and rs60/rs70/rs91 (GENE2), everything else is ultra-rare. Consequences:
#   - MAC25..MAC500 are all-ones on both genes, i.e. duplicates of unw;
#   - MAC1/3/5/10 select markers that are ALL collapsed at cutoff 10, so those
#     sets test the pseudo-marker alone (design doc, "What enters the test");
#   - with --MACCutoff_to_CollapseUltraRare=0 the whole MAC ladder is live on
#     real markers, which is the paper's unbinned quantitative run.
# 0.001 (MAC <= 2) would leave NO real marker at cutoff 10, which is why the
# tests and goldens pin --flexRV_maxMAF=0.01.

FX_G <- c("--bedFile=./input/genotype_100markers.bed",
          "--bimFile=./input/genotype_100markers.bim",
          "--famFile=./input/genotype_100markers.fam",
          "--AlleleOrder=alt-first")
FX_SG <- c(
  paste0("--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
         "randomMarkersUsed.sparseGRM.mtx"),
  paste0("--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_",
         "1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
)
FX_BIN <- c("--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
            "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt")
# the quantitative model runs WITHOUT the sparse GRM: its variance-ratio file
# has only `null` rows. A function, not a value: helper files are sourced in
# name order and helper-saige.R (which defines saige_goldens) comes after this
# one.
fx_quant <- function() {
  gold <- get("saige_goldens", envir = parent.frame())
  c(paste0("--GMMATmodelFile=", file.path(gold, "step1", "quant_sparseGRM.rda")),
    paste0("--varianceRatioFile=", file.path(gold, "step1", "quant_sparseGRM.varianceRatio.txt")))
}
FX_COMMON <- c("--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
               "--groupFile=./input/group_new_snpid_flexrv.txt")
# the FlexRV run proper: one annotation group, one max MAF, burden
FX_FLEX <- c("--flexRV_score=REVEL", "--flexRV_maxMAF=0.01",
             "--annotation_in_groupTest=lof:missense")
FX_BURDEN <- c("--r.corr=1", "--is_fastTest=FALSE")

# helper-saige.R's functions, bound here at top level: lintr cannot see them
# from inside a function body (the same aliasing test-step2-region-threads.R
# uses). This file is named so that it sorts AFTER helper-saige.R: testthat
# sources helpers in name order, and a helper-flexrv.R failed to load.
fx_skip <- skip_unless_level
fx_cli <- run_saige_cli
fx_read <- read_saige_table
fx_cli_ok <- expect_cli_ok

# One CLI run per distinct argument vector, cached per worker process.
fx_cache <- new.env(parent = emptyenv())
fx_run <- function(args, tag = "fx") {
  fx_skip("standard")
  key <- paste(args, collapse = " ")
  if (!is.null(fx_cache[[key]])) return(fx_cache[[key]])
  out <- file.path(tempfile(paste0("step2_", tag, "_")), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  log <- tempfile(fileext = ".log")
  # system2 goes through a shell: a transform spec's ';' would end the command
  args <- ifelse(grepl(";", args, fixed = TRUE), shQuote(args), args)
  status <- fx_cli("step2_SPAtests.R", c(args, paste0("--SAIGEOutputFile=", out)), log = log)
  region <- NULL
  if (identical(as.integer(status), 0L) && file.exists(out)) {
    region <- fx_read(out)
  }
  res <- list(status = status, log = readLines(log, warn = FALSE), path = out, region = region)
  fx_cache[[key]] <- res
  res
}
fx_ok <- function(res) {
  fx_cli_ok(res$status)
}

# The default 192-set binary BURDEN run (golden P58's command).
fx_burden_bin <- function() {
  fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN), "bin")
}
# The paper's quantitative run: no ultra-rare collapse (golden P59's command).
fx_burden_quant <- function() {
  fx_run(c(FX_G, fx_quant(), FX_COMMON, FX_FLEX, FX_BURDEN, "--MACCutoff_to_CollapseUltraRare=0"),
         "quant")
}

# Numeric columns of two region tables compared cell for cell at a relative
# tolerance, rows matched by (Region, Group, max_MAF, key) where `key` names
# the row within the cell (a Weight label, or NULL for the legacy schema).
fx_num_cols <- c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden", "SE_Burden")
fx_expect_rows_equal <- function(a, b, cols = NULL, tol = 1e-10, label = "rows") {
  if (is.null(cols)) {
    cols <- intersect(fx_num_cols, intersect(colnames(a), colnames(b)))
  }
  testthat::expect_identical(nrow(a), nrow(b), label = paste(label, "row count"))
  for (cc in cols) {
    x <- suppressWarnings(as.numeric(a[[cc]]))
    y <- suppressWarnings(as.numeric(b[[cc]]))
    testthat::expect_identical(is.na(x), is.na(y), label = paste(label, cc, "NA pattern"))
    testthat::expect_equal(x[!is.na(x)], y[!is.na(y)], tolerance = tol, label = paste(label, cc))
  }
}
