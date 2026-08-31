# Shared scaffolding for the --nThreads marker-scan tests, split across
# test-step2-threads.R (partition/merge units, P1 + P4, and the refusal
# guards) and test-step2-threads-goldens.R (P7 + P11 and the repeat-identity
# check). In a helper so both files can land in different parallel workers.
#
# Nothing runs at source time; th_run() is only called from inside test_that
# blocks that gate on skip_unless_level("standard").

# ---- end-to-end: threaded runs reproduce the serial goldens -----------------

TH_MODEL_BIN <- function() c(
  paste0("--GMMATmodelFile=",
         file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
  paste0("--varianceRatioFile=",
         file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt"))
)
TH_MODEL_QUANT <- function() c(
  paste0("--GMMATmodelFile=",
         file.path(saige_goldens, "step1", "quant_sparseGRM.rda")),
  paste0("--varianceRatioFile=",
         file.path(saige_goldens, "step1", "quant_sparseGRM.varianceRatio.txt"))
)
TH_PLINK_G <- c("--bedFile=./input/genotype_100markers.bed",
                "--bimFile=./input/genotype_100markers.bim",
                "--famFile=./input/genotype_100markers.fam",
                "--AlleleOrder=alt-first")
TH_PLINK_22CHR <- "./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr"

# The four marker goldens under test, each with make_goldens.sh's exact
# arguments minus --SAIGEOutputFile:
#   P1  ER path (minMAC=1, 70 markers)
#   P4  conditional, quantitative (COW-shared conditioning factors)
#   P7  multi-chunk, 5858 markers (chunking happens inside each shard)
#   P11 the AF>0.5 path, post-Stage-B
th_workloads <- list(
  P1 = list(golden = "plink_marker_binary.txt", args = function() c(
    TH_PLINK_G, "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    TH_MODEL_BIN(), "--is_output_moreDetails=TRUE")),
  P4 = list(golden = "plink_marker_cond_quantitative.txt", args = function() c(
    TH_PLINK_G, "--chrom=1", "--minMAF=0", "--minMAC=20", "--LOCO=FALSE",
    TH_MODEL_QUANT(), "--is_output_moreDetails=TRUE", "--condition=rs13,rs79")),
  P7 = list(golden = "plink_marker_chunked.txt", args = function() c(
    paste0("--bedFile=", TH_PLINK_22CHR, ".bed"),
    paste0("--bimFile=", TH_PLINK_22CHR, ".bim"),
    paste0("--famFile=", TH_PLINK_22CHR, ".fam"),
    "--AlleleOrder=alt-first",
    "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    TH_MODEL_BIN(), "--markers_per_chunk=1000")),
  P11 = list(golden = "plink_marker_reffirst.txt", args = function() c(
    "--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=ref-first",
    "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    TH_MODEL_BIN(), "--is_output_moreDetails=TRUE"))
)

th_run <- function(key, nThreads, tag = "") {
  out <- file.path(tempfile(paste0("step2_th_", key, "_", nThreads, tag, "_")),
                   "marker.txt")
  dir.create(dirname(out), recursive = TRUE)
  status <- run_saige_cli("step2_SPAtests.R",
                          c(th_workloads[[key]]$args(),
                            paste0("--SAIGEOutputFile=", out),
                            paste0("--nThreads=", nThreads)))
  expect_cli_ok(status)
  out
}

# nThreads = 64 exceeds the 70-marker P1/P11 sets only mildly, but exceeds the
# marker count of no workload here; it is included to exercise the clamp on
# P1/P11 (70 markers, 64 shards) and heavy oversubscription elsewhere.
TH_NTHREADS <- c(1, 2, 3, 4, 7, 64)
