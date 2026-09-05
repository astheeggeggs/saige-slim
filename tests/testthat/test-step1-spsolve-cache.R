# LEDGER #79 FIXED: the sparse-GRM direct solve keeps ONE SuperLU factorisation
# of Sigma = tau0/W + tau1*GRM alive while (W, tau) are unchanged, instead of
# rebuilding and refactorising the identical matrix for every right-hand side.
# Sigma depends only on (W, tau), which are constant across the right-hand
# sides of one getCoefficients call (the working vector, then each covariate
# column) and across the trace estimator's random vectors, so the number of
# factorisations a fit needs is the number of distinct (W, tau) pairs it
# visits, not the number of solves.
#
# The arithmetic is unchanged -- spsolve() at the default options and
# spsolve_factoriser both run get_perm_c, sp_preorder, gstrf, gstrs -- which is
# why every step-1 golden stays byte-identical (checked by the full tier and
# by hand when this landed). What this file gates is the SAVING: step 1 now
# reports "Sparse Sigma direct solves: S solves with F factorisations", and F
# must be a small fraction of S. At unpatched HEAD the line does not exist.
#
# Tier: standard (one sparse-GRM step-1 fit, the binary_sparseGRM_a golden's
# command, ~1 min).

test_that("LEDGER #79 FIXED: a sparse-GRM fit factorises Sigma far fewer times than it solves", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  prefix <- file.path(tempfile("step1_spsolve_"), "fit")
  dir.create(dirname(prefix), recursive = TRUE)
  st <- run_saige_cli("step1_fitNULLGLMM.R", c(
    "--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
    paste0("--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed",
           ".sparseGRM.mtx"),
    paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
           "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
    "--useSparseGRMtoFitNULL=TRUE",
    "--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",
    "--phenoCol=y_binary", "--covarColList=x1,x2,a9,a10", "--qCovarColList=a9,a10",
    "--sampleIDColinphenoFile=IID", "--traitType=binary",
    "--IsOverwriteVarianceRatioFile=TRUE", "--isCateVarianceRatio=TRUE",
    paste0("--outputPrefix=", prefix)
  ))
  expect_cli_ok(st)
  log <- readLines(attr(st, "log"), warn = FALSE)
  line <- grep("^Sparse Sigma direct solves: ", log, value = TRUE)
  expect_length(line, 1L)
  pat <- "^Sparse Sigma direct solves: ([0-9]+) solves with ([0-9]+) factorisations$"
  m <- regmatches(line, regexec(pat, line))[[1]]
  expect_length(m, 3L)
  solves <- as.integer(m[2])
  facts <- as.integer(m[3])
  # every distinct (W, tau) is factorised once; the fit visits one per AI-REML
  # iteration plus a few for the initial tau, and solves the working vector,
  # five covariate columns and the trace probes at each
  expect_gt(solves, 0L)
  expect_gt(facts, 0L)
  expect_lt(facts * 4L, solves)
  # The fitted model is the committed golden's, compared through the step-1
  # helpers rather than a hand-rolled bound: they are byte-exact on the golden
  # platform and fall back to the harness tolerances off it (TOL_THETA 1e-2 for
  # theta and the variance ratios, TOL_OTHER 1e-3 for the other fitted fields).
  #
  # This block used to assert rel_err < 1e-6 and byte-identity of the
  # variance-ratio file against a golden generated on another machine. That is a
  # cross-platform reproduction claim the estimator does not make: it passed on
  # macOS and failed on Linux CI at rel_err 2.6e-6, four orders of magnitude
  # inside the harness tolerance. It also bypassed SAIGE_COMPARE_MODE, so
  # prepush.sh's two modes both went green while CI went red.
  # helper-step1multi.R:16-26 records the same lesson being learned once
  # already, at 1e-5.
  gold <- file.path(saige_goldens, "step1", "binary_sparseGRM_a")
  expect_step1_matches(paste0(gold, ".rda"), paste0(prefix, ".rda"))
  expect_step1_matches(paste0(gold, ".varianceRatio.txt"),
                       paste0(prefix, ".varianceRatio.txt"))
})
