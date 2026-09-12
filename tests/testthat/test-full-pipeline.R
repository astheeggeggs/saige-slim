# Full tier (SAIGE_TEST_LEVEL=full, same-platform only): regenerate everything
# through the harness scripts and compare against the committed goldens, plus
# the step-1 determinism double-run and the permuted-phenotype negative
# control. Thin shell-outs so the harness scripts stay the single source of
# truth; expect ~20+ minutes.

run_bash <- function(script, args) {
  log <- tempfile(fileext = ".log")
  status <- suppressWarnings(system2("bash", c(file.path(saige_repo, "harness", script), args),
                                     stdout = log, stderr = log))
  attr(status, "log") <- log
  status
}

test_that("full regeneration matches goldens (make_goldens + compare_goldens)", {
  skip_unless_level("full")
  skip_if_no_goldens()
  skip_if_not(compare_wants_byte(),
              "full tier is byte-exact and therefore same-platform only")
  cand <- tempfile("candidate_")
  expect_cli_ok(run_bash("make_goldens.sh", cand))
  expect_cli_ok(run_bash("compare_goldens.sh", c(saige_goldens, cand)))
  unlink(cand, recursive = TRUE)
})

test_that("step1 is bitwise deterministic under the pinned env", {
  skip_unless_level("full")
  outs <- character(2)
  for (i in 1:2) {
    prefix <- file.path(tempfile("det_"), "run")
    dir.create(dirname(prefix), recursive = TRUE)
    status <- run_saige_cli("step1_fitNULLGLMM.R", c(
      "--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
      "--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
      "--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt",
      "--useSparseGRMtoFitNULL=TRUE",
      "--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",
      "--phenoCol=y_binary", "--covarColList=x1,x2,a9,a10", "--qCovarColList=a9,a10",
      "--sampleIDColinphenoFile=IID", "--traitType=binary",
      "--IsOverwriteVarianceRatioFile=TRUE", "--isCateVarianceRatio=TRUE",
      paste0("--outputPrefix=", prefix)))
    expect_cli_ok(status)
    outs[i] <- prefix
  }
  expect_files_byte_identical(paste0(outs[1], ".rda"), paste0(outs[2], ".rda"))
  expect_files_byte_identical(paste0(outs[1], ".varianceRatio.txt"),
                              paste0(outs[2], ".varianceRatio.txt"))
})

test_that("permuted-phenotype negative control passes and matches golden bytes", {
  skip_unless_level("full")
  skip_if_no_goldens()
  cand <- tempfile("control_")
  # make_controls.sh runs check_null.R itself and exits nonzero on a
  # non-uniform null, so expect_cli_ok covers the uniformity check.
  expect_cli_ok(run_bash("make_controls.sh", cand))
  if (compare_wants_byte()) {
    for (f in c("pheno_permuted.txt", "null_binary.rda",
                "null_binary.varianceRatio.txt", "null_marker.txt")) {
      expect_files_byte_identical(file.path(saige_goldens, "control", f),
                                  file.path(cand, "control", f))
    }
  }
  unlink(cand, recursive = TRUE)
})
