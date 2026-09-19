# Regression guard for LEDGER #158 ON THE FINISHED BUILD (2026-09-18).
#
# Why this is NOT harness/skato_pvalue/nullrep/dump_arm.R: that script patches
# SKAT:::Met_SKAT_Get_Pvalue, which on the finished build is NO LONGER the funnel
# every SKAT-O cell passes -- get_SKAT_pvalue routes to skatoScreenKernel() when
# the screen is on (default 0.05), so a screened cell never reaches SKAT's kernel
# and would be invisible to that patch. Re-running it as-is would silently score a
# SMALLER, BIASED cell set (the unscreened tail only) and read as a pass.
#
# So patch SAIGE:::get_SKAT_pvalue instead -- the one funnel, screened or not
# (SAIGE_SPATest_Region.R's skat_fun default and skatoMemoKernel's kernel default
# both resolve it at call time, so the patch is seen by every caller) -- and
# record, per cell, BOTH:
#   p          the SHIPPED answer (screen 0.05 + saddlepoint tail below 1e-5)
#   p_upstream SKAT's own kernel on the SAME (Score, Phi, r.corr)
# A screen MISS is then a cell with p > threshold and p_upstream <= threshold,
# countable directly rather than inferred. nThreads = 1 so canParallelSkato() is
# FALSE and no fork carries the cells away from this process.
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]; rds <- args[2]; D <- args[3]; M <- args[4]
suppressPackageStartupMessages({ library(SAIGE); library(SKAT) })

orig <- get("get_SKAT_pvalue", asNamespace("SAIGE"))
skatNs <- asNamespace("SKAT")
metOpt <- get("Met_SKAT_Get_Pvalue", skatNs)
cells <- vector("list", 20000); k <- 0L

patched <- function(Score, Phi, r.corr, regionTestType) {
  t0 <- proc.time()[["elapsed"]]
  shipped <- orig(Score, Phi, r.corr, regionTestType)
  t1 <- proc.time()[["elapsed"]]
  up <- try(metOpt(Score = Score, Phi = Phi, r.corr = r.corr,
                   method = "optimal.adj", Score.Resampling = NULL), silent = TRUE)
  t2 <- proc.time()[["elapsed"]]
  bad <- inherits(up, "try-error")
  k <<- k + 1L
  cells[[k]] <<- list(
    Score = Score, Phi = Phi, r.corr = r.corr,
    p = shipped$Pvalue_SKATO,
    p_skato_shipped = shipped$Pvalue_SKATO,
    p_skat_col = shipped$Pvalue_SKAT,
    p_burden_col = shipped$Pvalue_Burden,
    se_burden = shipped$SE_Burden,
    method = if (is.null(shipped$skato_method)) NA_character_ else shipped$skato_method,
    screen_p = if (is.null(shipped$skato_screen_p)) NA_real_ else shipped$skato_screen_p,
    p_upstream = if (bad) NA_real_ else up$p.value,
    p_each = if (bad) NULL else tryCatch(up$param$p.val.each, error = function(e) NULL),
    t = t1 - t0, t_upstream = t2 - t1)
  shipped
}
assignInNamespace("get_SKAT_pvalue", patched, ns = "SAIGE")

SPAGMMATtest(
  pgenFile = file.path(D, "chr20_sim_pgen.pgen"), pvarFile = file.path(D, "chr20_sim_pgen.pvar"),
  psamFile = file.path(D, "chr20_sim_pgen.psam"),
  chrom = "20", min_MAF = 0, min_MAC = 0.5,
  GMMATmodelFile = paste0(M, ".rda"), varianceRatioFile = paste0(M, ".varianceRatio.txt"),
  sparseGRMFile = file.path(D, "pedigree_truth.sparseGRM.mtx"),
  sparseGRMSampleIDFile = file.path(D, "pedigree_truth.sparseGRM.mtx.sampleIDs.txt"),
  LOCO = FALSE, is_Firth_beta = TRUE, pCutoffforFirth = 0.1, is_output_moreDetails = TRUE,
  is_fastTest = TRUE, is_output_markerList_in_groupTest = TRUE, is_single_in_groupTest = TRUE,
  maxMAF_in_groupTest = c(0.0001, 0.001, 0.01), groupFile = file.path(D, "groupFile.txt"),
  annotation_in_groupTest = c("high", "medium", "low", "high;medium", "high;medium;low"),
  nThreads = 1L, SAIGEOutputFile = out)

saveRDS(cells[seq_len(k)], rds)
cat("dumped", k, "cells\n")
