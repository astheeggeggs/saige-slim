# LEDGER #164 and #165 (2026-09-26): when the saddlepoint fails on a binary trait the
# fallback is the EXACT test, not the normal-approximation p, and every
# fallback is counted and named (R/SpaFallbacks.R). No bundled model can make
# a saddlepoint fail (min fitted probability 0.005, so LEDGER #163's overflow
# is out of reach), so the failure is PLANTED: SAIGE_SPA_PLANT_FAIL=1 declares
# every attempted saddlepoint non-converged (src/SAIGE_test.cpp). Two
# controls and one reconciliation:
#   negative -- without the plant the run prints "SPA fallbacks: 0", writes no
#               sidecar, and its output is the golden (byte, in byte mode);
#   positive -- with the plant every marker whose saddlepoint was attempted
#               (Is.SPA true in the golden and MAC above --max_MAC_for_ER) is
#               counted, listed, and its p is NOT the normal p;
#   two ways -- each planted marker's p equals the p of the SAME marker in a
#               run where --max_MAC_for_ER sends every marker to the exact
#               test natively. Same kernel, same per-marker seed, so the two
#               strings are identical, which is what makes "fell back to the
#               exact test" a checkable claim and not a label.
# The region path is covered on the FlexRV golden P58 (its MAC1 sets' ultra-
# rare pseudo-markers are the rows whose saddlepoint fires).

skip_unless_level("standard")
skip_if_no_goldens()

PLANT <- c(SAIGE_SPA_PLANT_FAIL = "1")

marker_args <- function(out, extra = character()) {
  c("--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first",
    paste0("--SAIGEOutputFile=", out),
    "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    paste0("--GMMATmodelFile=", file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
    paste0("--varianceRatioFile=", file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt")),
    "--is_output_moreDetails=TRUE", extra)
}

read_log <- function(status) readLines(attr(status, "log"), warn = FALSE)
fallback_line <- function(log) grep("^SPA fallbacks: ", log, value = TRUE)
fallback_count <- function(log) {
  as.integer(sub("^SPA fallbacks: ([0-9]+) .*$", "\\1", fallback_line(log)))
}

test_that("#164 negative control: a clean marker run reports 0 fallbacks and no sidecar", {
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", marker_args(out))
  expect_cli_ok(st)
  log <- read_log(st)
  expect_length(fallback_line(log), 1L)
  expect_identical(fallback_count(log), 0L)
  expect_false(file.exists(paste0(out, ".spaFallbacks.txt")))
  expect_output_matches_golden(file.path(saige_goldens, "step2", "plink_marker_binary.txt"), out)
})

test_that("#164 planted: every attempted saddlepoint falls back to the exact test, counted and named", {
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", marker_args(out), env = PLANT)
  expect_cli_ok(st)
  log <- read_log(st)
  gold <- read_saige_table(file.path(saige_goldens, "step2", "plink_marker_binary.txt"))
  res <- read_saige_table(out)
  # which golden rows had a saddlepoint attempted: Is.SPA true AND not the ER
  # route (MAC above the default --max_MAC_for_ER of 4; LEDGER #80 makes ER
  # markers read Is.SPA true too)
  # read_saige_table keeps every column as printed (character)
  ac <- as.numeric(gold$AC_Allele2)
  n2 <- 2 * (as.numeric(gold$N_case) + as.numeric(gold$N_ctrl))
  mac <- pmin(ac, n2 - ac)
  spa_rows <- which(tolower(gold$Is.SPA) == "true" & mac > 4)
  expect_gt(length(spa_rows), 0)   # P1 has two such markers; the region test below has 16
  n <- fallback_count(log)
  expect_identical(n, length(spa_rows))
  side <- paste0(out, ".spaFallbacks.txt")
  expect_true(file.exists(side))
  tab <- read_saige_table(side)
  expect_identical(nrow(tab), length(spa_rows))
  expect_setequal(tab$MarkerID, gold$MarkerID[spa_rows])
  expect_true(all(tab$Kind == "marker"))
  # the normal p was recorded as NOT reported, and is not what the file carries
  expect_identical(as.character(tab$p_normal_not_reported),
                   as.character(res$p.value.NA[match(tab$MarkerID, res$MarkerID)]))
  expect_identical(as.character(tab$p_exact_reported),
                   as.character(res$p.value[match(tab$MarkerID, res$MarkerID)]))
  expect_false(any(res$p.value[spa_rows] == res$p.value.NA[spa_rows]))
  # rows whose saddlepoint was never attempted are untouched
  other <- setdiff(seq_len(nrow(gold)), spa_rows)
  expect_identical(as.character(res$p.value[other]), as.character(gold$p.value[other]))
})

test_that("#164 two ways: the planted fallback p equals the native exact-test p", {
  out_p <- tempfile(fileext = ".txt")
  out_e <- tempfile(fileext = ".txt")
  st_p <- run_saige_cli("step2_SPAtests.R", marker_args(out_p), env = PLANT)
  # --max_MAC_for_ER=200 is the value golden P33 proved on this fixture, so the
  # native exact test covers every planted marker with MAC <= 200
  st_e <- run_saige_cli("step2_SPAtests.R", marker_args(out_e, "--max_MAC_for_ER=200"))
  expect_cli_ok(st_p); expect_cli_ok(st_e)
  # the native ER run itself must not fall back (an ER marker attempts no saddlepoint)
  expect_identical(fallback_count(read_log(st_e)), 0L)
  tab <- read_saige_table(paste0(out_p, ".spaFallbacks.txt"))
  rp <- read_saige_table(out_p); re <- read_saige_table(out_e)
  ip <- match(tab$MarkerID, rp$MarkerID); ie <- match(tab$MarkerID, re$MarkerID)
  expect_false(anyNA(c(ip, ie)))
  # The fallback runs the exact test only in its enumeration regime (<= 13
  # carriers); above that it reports NA and the sidecar's Note says so. The
  # two P1 markers whose saddlepoint fires are both above it, so on this
  # fixture the marker path exercises the NA branch and the region path (the
  # cross-arm test below, 4 carriers) exercises the identity with native ER.
  ceiling_rows <- grepl("above the fallback's exact-test ceiling", tab$Note)
  expect_true(all(rp$p.value[ip][ceiling_rows] == "NA"))
  expect_true(all(grepl("^spa: planted", tab$Note)))
  ran <- !ceiling_rows
  if (any(ran)) {
    expect_identical(as.character(rp$p.value[ip][ran]), as.character(re$p.value[ie][ran]))
  }
})

test_that("#164 region path: the planted fallback reaches the FlexRV ultra-rare rows and the BURDEN cells", {
  out <- tempfile(fileext = ".txt")
  args <- c("--bedFile=./input/genotype_100markers.bed",
            "--bimFile=./input/genotype_100markers.bim",
            "--famFile=./input/genotype_100markers.fam",
            "--AlleleOrder=alt-first",
            paste0("--SAIGEOutputFile=", out),
            "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
            "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
            "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
            "--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
            "--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt",
            "--groupFile=./input/group_new_snpid_flexrv.txt",
            "--flexRV_score=REVEL", "--flexRV_maxMAF=0.01",
            "--annotation_in_groupTest=lof:missense",
            "--r.corr=1", "--is_fastTest=FALSE", "--is_single_in_groupTest=TRUE",
            # LEDGER #166 made every set ER-eligible at MAC <= 4, so the MAC1
            # sets no longer reach the saddlepoint natively; with the cutoff at
            # 0 they do, and the plant then exercises the fallback as before
            "--max_MAC_for_ER=0")
  st0 <- run_saige_cli("step2_SPAtests.R", args)
  expect_cli_ok(st0)
  expect_identical(fallback_count(read_log(st0)), 0L)
  # (not the golden: P58 runs at the default cutoff, where the MAC1 sets take the exact test)
  st <- run_saige_cli("step2_SPAtests.R", args, env = PLANT)
  expect_cli_ok(st)
  log <- read_log(st)
  side <- paste0(out, ".spaFallbacks.txt")
  expect_true(file.exists(side))
  tab <- read_saige_table(side)
  expect_identical(fallback_count(log), nrow(tab))
  expect_setequal(unique(tab$Kind), c("ultra_rare", "burden"))
  # every recorded ultra-rare pseudo-marker is a row of the sidecar whose
  # saddlepoint fired in the clean run (the MAC1 sets), and its cell labels
  # name a row of the results file
  single <- read_saige_table(paste0(out, ".singleAssoc.txt"))
  ur <- tab[tab$Kind == "ultra_rare", ]
  expect_true(all(ur$MarkerID %in% single$MarkerID))
  expect_true(all(grepl("_MAC1$", ur$Weight)))
  res <- read_saige_table(out)
  bu <- tab[tab$Kind == "burden", ]
  key_tab <- paste(bu$Region, bu$Group, bu$Weight)
  key_res <- paste(res$Region, res$Group, res$Weight)
  expect_true(all(key_tab %in% key_res))
  # the reported burden p on those cells is the exact p, not the normal p
  m <- match(key_tab, key_res)
  expect_identical(as.character(bu$p_exact_reported), as.character(res$Pvalue_Burden[m]))
  expect_false(any(as.character(bu$p_normal_not_reported) == as.character(res$Pvalue_Burden[m])))
})

test_that("#164 cross-arm: the planted FlexRV set agrees with the SKAT-O arm's native exact test", {
  # The AoU shape end to end. The FlexRV arm's `unw_MAC1` set collapses the
  # MAC <= 1 markers with unit weights; the SKAT-O arm's beta set at max MAF
  # 0.0005 (MAC <= 1 on 1,000 samples) collapses the SAME markers, and being
  # a beta set it is ER-eligible, so at MAC 4 it takes the exact test
  # natively. With the plant, the FlexRV set's saddlepoint fails and falls
  # back to that exact test: the two arms must now print the same p for the
  # same vector, which is exactly what they did NOT do before this fix.
  base <- c("--bedFile=./input/genotype_100markers.bed",
            "--bimFile=./input/genotype_100markers.bim",
            "--famFile=./input/genotype_100markers.fam",
            "--AlleleOrder=alt-first",
            "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
            "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
            "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
            "--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
            "--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt",
            "--groupFile=./input/group_new_snpid_flexrv.txt",
            "--annotation_in_groupTest=lof:missense",
            "--is_fastTest=FALSE", "--is_single_in_groupTest=TRUE")
  out_f <- tempfile(fileext = ".txt")
  st_f <- run_saige_cli("step2_SPAtests.R",
                        c(base, paste0("--SAIGEOutputFile=", out_f), "--flexRV_score=REVEL",
                          "--flexRV_maxMAF=0.01", "--r.corr=1", "--max_MAC_for_ER=0"), env = PLANT)
  # LEDGER #166: at the DEFAULT cutoff and with no plant, the same set takes the
  # exact test natively, like the SKAT-O arm's beta set -- the AoU discrepancy
  # closed without any fallback
  out_n <- tempfile(fileext = ".txt")
  st_n <- run_saige_cli("step2_SPAtests.R",
                        c(base, paste0("--SAIGEOutputFile=", out_n), "--flexRV_score=REVEL",
                          "--flexRV_maxMAF=0.01", "--r.corr=1"))
  out_s <- tempfile(fileext = ".txt")
  st_s <- run_saige_cli("step2_SPAtests.R",
                        c(base, paste0("--SAIGEOutputFile=", out_s), "--maxMAF_in_groupTest=0.0005"))
  expect_cli_ok(st_f); expect_cli_ok(st_s); expect_cli_ok(st_n)
  expect_identical(fallback_count(read_log(st_s)), 0L)
  expect_identical(fallback_count(read_log(st_n)), 0L)
  sn <- read_saige_table(paste0(out_n, ".singleAssoc.txt"))
  rn <- sn[sn$MarkerID == "GENE1:lof;missense:0.01:REVEL_unw_MAC1", ]
  expect_identical(nrow(rn), 1L)
  sf <- read_saige_table(paste0(out_f, ".singleAssoc.txt"))
  ss <- read_saige_table(paste0(out_s, ".singleAssoc.txt"))
  rf <- sf[sf$MarkerID == "GENE1:lof;missense:0.01:REVEL_unw_MAC1", ]
  rs <- ss[grepl("^GENE1:lof;missense:0\\.0005", ss$MarkerID), ]
  expect_identical(nrow(rf), 1L); expect_identical(nrow(rs), 1L)
  # same vector: same score statistic and same normal p
  expect_identical(rf$Tstat, rs$Tstat)
  expect_identical(rf$p.value.NA, rs$p.value.NA)
  # the SKAT-O arm's p is the native exact test; the planted FlexRV arm now matches it
  expect_true(tolower(rs$Is.SPA) == "true")
  expect_identical(rf$p.value, rs$p.value)
  expect_identical(rn$p.value, rs$p.value)        # native (#166) == SKAT-O arm's exact p
  expect_true(tolower(rn$Is.SPA) == "true")
  side <- read_saige_table(paste0(out_f, ".spaFallbacks.txt"))
  expect_true("GENE1:lof;missense:0.01:REVEL_unw_MAC1" %in% side$MarkerID)
})

test_that("#165: --max_MAC_for_ER above the verified ceiling is refused", {
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", marker_args(out, "--max_MAC_for_ER=201"))
  expect_false(st == 0)
  expect_true(any(grepl("max_MAC_use_ER should be a numeric value less than or equal to 200",
                        read_log(st))))
})
