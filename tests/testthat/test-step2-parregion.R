# LEDGER #38 FIXED: a --X_PARregion token that is not "start-end" is refused
# by name instead of being recycled by rbind into a boundary row (a lone
# "5" became c(5, 5)) and handed to C++ silently. Tier: standard (one refusal
# that fails before any marker is read; the well-formed default is exercised
# by test-step2-xnonpar-chrom.R).

test_that("a malformed PAR token is refused by name", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  m <- file.path(saige_goldens, "step1", "binary_sparseGRM_a")
  e <- new.env(); load(paste0(m, ".rda"), envir = e)
  males <- tempfile("males_", fileext = ".txt"); writeLines(as.character(e$modglmm$sampleID)[1:10], males)
  for (bad in c("60001,154931044-155260560", "60001-abc", "700-600")) {
    st <- run_saige_cli("step2_SPAtests.R", c(
      "--bedFile=./input/genotype_100markers.bed", "--bimFile=./input/genotype_100markers.bim",
      "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
      paste0("--SAIGEOutputFile=", tempfile(fileext = ".txt")), "--chrom=1", "--minMAF=0", "--minMAC=1",
      "--LOCO=FALSE", paste0("--GMMATmodelFile=", m, ".rda"), paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
      "--is_rewrite_XnonPAR_forMales=TRUE", paste0("--sampleFile_male=", males), paste0("--X_PARregion=", bad)))
    expect_false(identical(as.integer(st), 0L))
    lg <- readLines(attr(st, "log"), warn = FALSE)
    expect_true(any(grepl("X_PARregion token '.*' is not of the form start-end", lg)))
  }
})
