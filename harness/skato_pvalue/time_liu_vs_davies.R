# Times SKAT's Davies-based optimal.adj against modified Liu (liu.mod) on every
# SKAT-O cell of a sim50k run; the Davies result is what the run keeps, so the
# output file is the unpatched one. One record per kernel call.
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]; rec <- args[2]; D <- args[3]
suppressPackageStartupMessages(library(SAIGE))
orig <- get("Met_SKAT_Get_Pvalue", envir = asNamespace("SKAT"))
con <- file(rec, open = "w")
writeLines("q\tt_davies\tt_liu\tp_davies\tp_liu\tp_skat_davies\tp_skat_liu\tp_burden_davies\tp_burden_liu\tliu_err", con)
timed <- function(Score, Phi, r.corr, method, Score.Resampling = NULL, ...) {
  t0 <- proc.time()[["elapsed"]]
  d <- orig(Score = Score, Phi = Phi, r.corr = r.corr, method = method, Score.Resampling = Score.Resampling, ...)
  t1 <- proc.time()[["elapsed"]]
  l <- try(orig(Score = Score, Phi = Phi, r.corr = r.corr, method = "liu.mod", Score.Resampling = Score.Resampling, ...), silent = TRUE)
  t2 <- proc.time()[["elapsed"]]
  err <- inherits(l, "try-error")
  g <- function(x, i) { v <- try(x$p.val.each[i], silent = TRUE); if (inherits(v, "try-error") || length(v) == 0) NA else v }
  nr <- length(r.corr)
  writeLines(sprintf("%d\t%.6f\t%.6f\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%d",
    length(Score), t1 - t0, t2 - t1, d$p.value, if (err) NA else l$p.value,
    g(d, 1), if (err) NA else g(l, 1), g(d, nr), if (err) NA else g(l, nr), as.integer(err)), con)
  d
}
assignInNamespace("Met_SKAT_Get_Pvalue", timed, ns = "SKAT")
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
close(con)
