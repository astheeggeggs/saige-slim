# Shared by test-step1-pheno-scale.R (LEDGER #65) and
# test-step1-vr-sparse-scale.R (LEDGER #97): the committed phenotype file with
# y_quantitative multiplied by `scale`, everything else untouched. The .bed,
# sparse GRM and covariates are the quant_sparseGRM golden's; invNormalize is
# deliberately NOT used by the tests that read this (it would undo the scaling
# and hide the bug).
ps_pheno_file <- function(scale) {
  src <- file.path(saige_extdata, "input",
                   "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt")
  d <- read.table(src, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  d$y_quantitative <- d$y_quantitative * scale
  out <- tempfile(sprintf("pheno_x%g_", scale), fileext = ".txt")
  write.table(d, out, quote = FALSE, row.names = FALSE, sep = " ")
  out
}
