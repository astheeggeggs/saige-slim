# Regression test for LEDGER.md #16 (FIXED in Phase 3e): a mid-stream request
# for bim row 0 used to perform no seek and decode whatever record the cursor
# sat on (conditioning on rs80 made monomorphic rs1 appear carrying rs81's
# AC=314). The reader now seeks unconditionally: rs1 must be absent from BOTH
# runs (it is truly monomorphic and filtered).
#
# The conditioning marker is rs79 now, not rs80: rs80 is itself monomorphic in
# the analysed samples and LEDGER #73 (FIXED 2026-09-05) refuses such a
# conditioning marker by name. Any marker read before the loop puts the cursor
# mid-stream, which is all the original construction needed.

test_that("LEDGER #16 fixed: mid-stream bim-row-0 reads return the true record", {
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
                              c(base_args, "--condition=rs79",
                                paste0("--SAIGEOutputFile=", out2))))
  cond <- read_saige_table(out2)
  expect_false("rs1" %in% cond$MarkerID)  # no phantom: true (monomorphic) rs1 is filtered
  # and the conditioning did not disturb any other marker's unconditional stats
  shared <- intersect(colnames(fresh), colnames(cond))
  expect_identical(fresh[fresh$MarkerID == "rs79", shared],
                   cond[cond$MarkerID == "rs79", shared])
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
