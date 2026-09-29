# LEDGER #173 (DECIDED 2026-09-29: upstream's arithmetic stays) and its safety
# net. On the marker path under the default --is_noadjCov=TRUE the score is the
# CENTRED genotype's and the saddlepoint's CGF the covariate-PROJECTED one's,
# so the fast kernel's Gaussian remainder, var2(centred) minus the projected
# carriers' part, is not a variance (test-step2-spa-dump.R characterizes it).
# It could in principle go negative; then the dense kernel, which needs no
# remainder, runs instead and the run counts it. It has never been observed,
# so the net is exercised by a plant: SAIGE_SPA_PLANT_NEG_REMAINDER=1 declares
# every fast-kernel remainder bad.
#   negative control  P1 as shipped: the golden, and no `SPA remainder` line;
#   positive control  P1 planted: the line counts the saddlepoints that would
#                     have taken the fast kernel, every dump names the dense
#                     kernel and replays to the printed p, the rows that move
#                     are saddlepoint rows and no others, and the two kernels
#                     agree on the tail to within 25 %.
skip_if_no_goldens()

remainder_args <- function(out) {
  c("--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first",
    paste0("--SAIGEOutputFile=", out),
    "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    paste0("--GMMATmodelFile=", file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
    paste0("--varianceRatioFile=",
           file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt")),
    "--is_output_moreDetails=TRUE")
}
golden_marker <- file.path(saige_goldens, "step2", "plink_marker_binary.txt")
remainder_line <- function(log) grep("^SPA remainder: ", log, value = TRUE)

test_that("as shipped no remainder is refused: the golden, and the run says nothing", {
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", remainder_args(out))
  expect_cli_ok(st)
  expect_output_matches_golden(golden_marker, out)
  expect_length(remainder_line(readLines(attr(st, "log"), warn = FALSE)), 0L)
})

test_that("planted, every fast-kernel saddlepoint runs the dense kernel and is counted", {
  out <- tempfile(fileext = ".txt")
  dir <- tempfile("spadump_")
  dir.create(dir)
  st <- run_saige_cli("step2_SPAtests.R", remainder_args(out),
                      env = c(SAIGE_SPA_PLANT_NEG_REMAINDER = "1", SAIGE_SPA_DUMP_DIR = dir,
                              SAIGE_SPA_DUMP_BELOW = "1"))
  expect_cli_ok(st)
  log <- readLines(attr(st, "log"), warn = FALSE)
  gold <- read_saige_table(golden_marker)
  new <- read_saige_table(out)
  expect_identical(dim(new), dim(gold))
  moved <- which(new$p.value != gold$p.value)
  expect_gt(length(moved), 0)
  n_line <- as.integer(sub("^SPA remainder: ([0-9]+) test.*$", "\\1", remainder_line(log)))
  expect_identical(n_line, length(moved))
  # only the saddlepoint's rows move, and only in what the saddlepoint sets
  expect_true(all(tolower(gold$Is.SPA[moved]) == "true"))
  same_cols <- setdiff(names(gold), c("p.value", "SE"))
  expect_identical(new[, same_cols], gold[, same_cols])
  expect_identical(new[-moved, ], gold[-moved, ])
  ratio <- as.numeric(new$p.value[moved]) / as.numeric(gold$p.value[moved])
  expect_true(all(abs(ratio - 1) < 0.25))
  # every converged saddlepoint is dumped, names the dense kernel and replays
  # to the p the run printed. P1 has two: one took the dense kernel already
  # (a common marker, more carriers than not) and is where it was; the other
  # is the one the plant moved.
  files <- list.files(dir, pattern = "^spaBelow_[0-9]+_[0-9]+\\.tsv$", full.names = TRUE)
  ac <- as.numeric(gold$AC_Allele2)
  n2 <- 2 * (as.numeric(gold$N_case) + as.numeric(gold$N_ctrl))
  spa_rows <- which(tolower(gold$Is.SPA) == "true" & pmin(ac, n2 - ac) > 4)
  expect_length(files, length(spa_rows))
  expect_true(all(moved %in% spa_rows))
  expect_lt(length(moved), length(spa_rows))
  printed <- as.numeric(new$p.value[spa_rows])
  for (f in files) {
    dm <- read_spa_dump(f)
    expect_identical(dm$h2[["kernel"]], "dense")
    rp <- replay_spa_dump(dm)
    expect_true(isTRUE(rp$Isconverge))
    expect_equal(rp$pvalue, spa_dump_num(dm$h2, "p"), tolerance = 1e-10)
    expect_true(any(abs(printed / rp$pvalue - 1) < 1e-5))
  }
})
