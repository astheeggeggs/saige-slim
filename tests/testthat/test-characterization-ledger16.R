# Characterization of LEDGER.md #16 (CONFIRMED bug, preserved bug-for-bug):
# a mid-stream request for bim row 0 performs no seek in the PLINK reader and
# reads whatever record the cursor sits on. Reachable via conditional analysis
# (conditioning markers are read before the marker loop).
#
# This test PINS THE WRONG BEHAVIOR on purpose: conditioning on rs80 parks the
# cursor on rs81's record (AC=314); the loop then requests rs1 (bim row 0,
# truly monomorphic) and reports it with rs81's genotypes. When the 3e fix
# lands, update this test together with the ledger entry: rs1 must then be
# absent from the conditional output too.

test_that("LEDGER #16: mid-stream bim-row-0 read returns the cursor's record", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  ids <- tempfile(); writeLines(c("rs1", "rs79", "rs80"), ids)
  qm <- c(
    paste0("--GMMATmodelFile=", file.path(saige_goldens, "step1", "quant_sparseGRM.rda")),
    paste0("--varianceRatioFile=", file.path(saige_goldens, "step1", "quant_sparseGRM.varianceRatio.txt")))
  base_args <- c(
    "--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first", "--chrom=1", "--minMAF=0", "--minMAC=0.5",
    "--LOCO=FALSE", qm, paste0("--idstoIncludeFile=", ids))

  out1 <- tempfile(fileext = ".txt")
  expect_cli_ok(run_saige_cli("step2_SPAtests.R",
                              c(base_args, paste0("--SAIGEOutputFile=", out1))))
  fresh <- read_saige_table(out1)
  expect_false("rs1" %in% fresh$MarkerID)  # fresh cursor: rs1 is monomorphic, filtered

  out2 <- tempfile(fileext = ".txt")
  expect_cli_ok(run_saige_cli("step2_SPAtests.R",
                              c(base_args, "--condition=rs80",
                                paste0("--SAIGEOutputFile=", out2))))
  cond <- read_saige_table(out2)
  rs1 <- cond[cond$MarkerID == "rs1", ]
  expect_identical(nrow(rs1), 1L)               # the phantom appears...
  expect_identical(rs1$AC_Allele2, "314")       # ...carrying rs81's allele count
})

test_that("LEDGER #10 fixed: marker run without null_noXadj VR rows stops cleanly", {
  skip_unless_level("standard")
  status <- run_saige_cli("step2_SPAtests.R", c(
    "--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first", "--chrom=1", "--minMAF=0", "--minMAC=1",
    "--LOCO=FALSE",
    "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
    "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
    paste0("--SAIGEOutputFile=", tempfile(fileext = ".txt"))))
  expect_false(identical(as.integer(status), 0L))
  lg <- readLines(attr(status, "log"), warn = FALSE)
  expect_true(any(grepl("no 'null_noXadj' rows", lg)))
})
