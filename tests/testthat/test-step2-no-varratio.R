# LEDGER #70 FIXED: an omitted --varianceRatioFile is now an R warning, not a
# bare cat(). The behaviour is unchanged on purpose -- every marker gets a
# variance ratio of exactly 1, which is what a step-1 run that skipped the
# estimation intends -- but the substitution scales var by 1/r, BETA by r and
# SE by sqrt(r) in EVERY row, so it has to be visible. This test pins both
# halves: the warning is raised, and the scaling is exactly what the row
# says, reconciled against golden P1 row by row using the ratios in the
# step-1 file. Tier: standard (one marker run).

test_that("no --varianceRatioFile: warns, and every row is P1 scaled by one of the file's ratios", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  m <- file.path(saige_goldens, "step1", "binary_sparseGRM_a")
  out <- tempfile(fileext = ".txt")
  args <- c("--bedFile=./input/genotype_100markers.bed", "--bimFile=./input/genotype_100markers.bim",
            "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
            paste0("--SAIGEOutputFile=", out), "--chrom=1", "--minMAF=0", "--minMAC=1",
            "--LOCO=FALSE", paste0("--GMMATmodelFile=", m, ".rda"),
            "--is_output_moreDetails=TRUE")
  st <- run_saige_cli("step2_SPAtests.R", args)
  expect_cli_ok(st)
  log <- readLines(attr(st, "log"), warn = FALSE)
  expect_true(any(grepl("varianceRatioFile is not specified, so a variance ratio of 1 is used", log, fixed = TRUE)))
  expect_true(any(grepl("^Warning message", log)))

  p1 <- read_saige_table(file.path(saige_goldens, "step2", "plink_marker_binary.txt"))
  d <- read_saige_table(out)
  expect_identical(d$MarkerID, p1$MarkerID)
  ratios <- utils::read.table(paste0(m, ".varianceRatio.txt"))$V1
  # var(P1) = r * var(here) for the ratio r of that marker's MAC category
  r_obs <- as.numeric(p1$var) / as.numeric(d$var)
  ok <- is.finite(r_obs)
  expect_gt(sum(ok), 60)
  nearest <- vapply(r_obs[ok], function(x) min(abs(x - ratios) / ratios), numeric(1))
  expect_true(all(nearest < 1e-5))
  expect_true(all(abs(r_obs[ok] - 1) > 0.01))                       # it is not a no-op
  # BETA scales by r and SE by sqrt(r), the other half of the row's claim
  b_ratio <- as.numeric(d$BETA) / as.numeric(p1$BETA)
  s_ratio <- (as.numeric(d$SE) / as.numeric(p1$SE))^2
  fin <- ok & is.finite(b_ratio) & is.finite(s_ratio) & p1$`Is.SPA` == "false" & as.numeric(p1$BETA) != 0
  expect_gt(sum(fin), 30)
  expect_true(all(abs(b_ratio[fin] - r_obs[fin]) < 1e-4))
  expect_true(all(abs(s_ratio[fin] - r_obs[fin]) < 1e-4))
})
