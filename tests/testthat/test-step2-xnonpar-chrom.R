# LEDGER #69 FIXED: --is_rewrite_XnonPAR_forMales doubled male genotypes on
# EVERY chromosome, because processMale_XnonPAR decides from the position
# alone and the chromosome was never consulted. The docstring told the user to
# feed it an X-only file instead of the code checking. It checks now.
#
# No bundled fixture has chromosome X, so the positive arm DERIVES one: the
# bundled .bim with its chromosome column rewritten to "X" (the .bed is indexed
# by row, so no genotype changes). Males are a derived list: every other model
# sample. Three runs, all otherwise P1's command on the step-1 golden model:
#   base : no flag, chromosome 1
#   auto : the flag, chromosome 1     -> must be byte-identical to base now
#   x    : the flag, chromosome "X"   -> male dosages doubled, must DIFFER
# Tier: standard.

xnp_env <- new.env()

xnp_setup <- function() {
  if (!is.null(xnp_env$male)) return(invisible(NULL))
  m <- file.path(saige_goldens, "step1", "binary_sparseGRM_a")
  e <- new.env(); load(paste0(m, ".rda"), envir = e)
  ids <- as.character(e$modglmm$sampleID)
  xnp_env$male <- tempfile("males_", fileext = ".txt")
  writeLines(ids[seq(1, length(ids), by = 2)], xnp_env$male)
  bim <- readLines(file.path(saige_extdata, "input", "genotype_100markers.bim"), warn = FALSE)
  xnp_env$bimX <- tempfile("chrX_", fileext = ".bim")
  writeLines(sub("^1\t", "X\t", bim), xnp_env$bimX)
  expect_true(all(grepl("^X\t", readLines(xnp_env$bimX))))
  invisible(NULL)
}

xnp_run <- function(flag, chrX) {
  xnp_setup()
  m <- file.path(saige_goldens, "step1", "binary_sparseGRM_a")
  out <- tempfile(fileext = ".txt")
  args <- c("--bedFile=./input/genotype_100markers.bed",
            paste0("--bimFile=", if (chrX) xnp_env$bimX else "./input/genotype_100markers.bim"),
            "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
            paste0("--SAIGEOutputFile=", out), paste0("--chrom=", if (chrX) "X" else "1"),
            "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
            paste0("--GMMATmodelFile=", m, ".rda"),
            paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
            "--is_output_moreDetails=TRUE")
  if (flag) args <- c(args, "--is_rewrite_XnonPAR_forMales=TRUE",
                      paste0("--sampleFile_male=", xnp_env$male))
  st <- run_saige_cli("step2_SPAtests.R", args)
  expect_cli_ok(st)
  list(out = out, log = readLines(attr(st, "log"), warn = FALSE))
}

test_that("the flag leaves an autosome alone, and doubles males on chromosome X", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  base <- xnp_run(FALSE, FALSE)
  auto <- xnp_run(TRUE, FALSE)
  expect_true(any(grepl("males are found in the test", auto$log, fixed = TRUE)))  # the flag was live
  expect_identical(readLines(auto$out), readLines(base$out))                       # ... and inert on chr 1

  x <- xnp_run(TRUE, TRUE)
  b <- read_saige_table(base$out); d <- read_saige_table(x$out)
  expect_identical(d$MarkerID, b$MarkerID)
  ac_b <- as.numeric(b$AC_Allele2); ac_x <- as.numeric(d$AC_Allele2)
  # a doubled male dosage raises the ALT count wherever a male carries it;
  # positions 1..100 are outside both PAR windows, so every marker qualifies
  expect_gt(sum(ac_x > ac_b), 50)
  expect_true(all(ac_x >= ac_b))
})
