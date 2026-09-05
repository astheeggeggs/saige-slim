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

# the first nPerChr markers of chr 1 and of chr 2 of the step-1 fixture, as a
# PLINK set under `out`
sm_two_chr_fixture <- function(out, nPerChr = 1000L) {
  src <- file.path(saige_extdata, "input",
                   "nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr")
  bim <- read.table(paste0(src, ".bim"), header = FALSE, stringsAsFactors = FALSE,
                    colClasses = "character")
  fam <- readLines(paste0(src, ".fam"))
  recBytes <- (length(fam) + 3L) %/% 4L
  keep <- c(which(bim[[1]] == "1")[seq_len(nPerChr)], which(bim[[1]] == "2")[seq_len(nPerChr)])
  stopifnot(!anyNA(keep))
  con <- file(paste0(src, ".bed"), "rb")
  bed <- readBin(con, "raw", 3L + nrow(bim) * recBytes)
  close(con)
  stopifnot(length(bed) == 3L + nrow(bim) * recBytes)
  sel <- unlist(lapply(keep - 1L, function(i) 3L + i * recBytes + seq_len(recBytes)))
  writeBin(c(bed[1:3], bed[sel]), paste0(out, ".bed"))
  write.table(bim[keep, ], paste0(out, ".bim"), quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  writeLines(fam, paste0(out, ".fam"))
  out
}
