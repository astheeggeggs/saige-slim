# Dump every SKAT-O kernel input (Score, Phi, r.corr) and SKAT's own answer to an RDS.
args <- commandArgs(trailingOnly = TRUE); out <- args[1]; rds <- args[2]; D <- args[3]
suppressPackageStartupMessages({ library(SAIGE); library(SKAT) })
orig <- get("Met_SKAT_Get_Pvalue", asNamespace("SKAT"))
cells <- vector("list", 5000); k <- 0L
patched <- function(Score, Phi, r.corr, method, Score.Resampling = NULL, ...) {
  t0 <- proc.time()[["elapsed"]]
  d <- orig(Score = Score, Phi = Phi, r.corr = r.corr, method = method, Score.Resampling = Score.Resampling, ...)
  t1 <- proc.time()[["elapsed"]]
  k <<- k + 1L; cells[[k]] <<- list(Score = Score, Phi = Phi, r.corr = r.corr, method = method, p = d$p.value, p_each = tryCatch(d$param$p.val.each, error = function(e) NULL), t = t1 - t0)
  d
}
assignInNamespace("Met_SKAT_Get_Pvalue", patched, ns = "SKAT")
SPAGMMATtest(
  pgenFile = file.path(D, "chr20_sim_pgen.pgen"), pvarFile = file.path(D, "chr20_sim_pgen.pvar"), psamFile = file.path(D, "chr20_sim_pgen.psam"),
  chrom = "20", min_MAF = 0, min_MAC = 0.5,
  GMMATmodelFile = file.path(D, "step1_y_binary.rda"), varianceRatioFile = file.path(D, "step1_y_binary.varianceRatio.txt"),
  sparseGRMFile = file.path(D, "pedigree_truth.sparseGRM.mtx"), sparseGRMSampleIDFile = file.path(D, "pedigree_truth.sparseGRM.mtx.sampleIDs.txt"),
  LOCO = FALSE, is_Firth_beta = TRUE, pCutoffforFirth = 0.1, is_output_moreDetails = TRUE, is_fastTest = TRUE,
  is_output_markerList_in_groupTest = TRUE, is_single_in_groupTest = TRUE,
  maxMAF_in_groupTest = c(0.0001, 0.001, 0.01), groupFile = file.path(D, "groupFile.txt"),
  annotation_in_groupTest = c("high", "medium", "low", "high;medium", "high;medium;low"),
  nThreads = 1L, SAIGEOutputFile = out)
saveRDS(cells[seq_len(k)], rds); cat("dumped", k, "cells\n")
