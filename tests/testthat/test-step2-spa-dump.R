# SAIGE_SPA_DUMP_BELOW (2026-09-28, the LEDGER #172 investigation). With
# SAIGE_SPA_DUMP_DIR set, a threshold p dumps every CONVERGED saddlepoint whose
# p is below it as <dir>/spaBelow_<pid>_<k>.tsv; the #168 failure dump cannot
# see a saddlepoint that converged to a wrong value, which is what three AoU
# FlexRV cells did (1e-53..1e-93 with Is.SPA true, exact answer ~1e-3). Every
# dump now carries a second header line naming the kernel that ran, both roots,
# the Gaussian remainder (NAmu, NAsigma) and the reported p, plus an `nz` column
# flagging the carriers, so a dump REPLAYS through the exported kernels, and a
# `dropped` column (LEDGER #174) flagging the pinned samples left out of the CGF
# (none on P1's model; test-step2-spa-pinned.R replays dumps that have them).
#   1. the switch cannot move a result: P1 with the dump on IS the golden;
#   2. at threshold 1 the dump count is the number of golden rows whose
#      saddlepoint converged (Is.SPA true, MAC above --max_MAC_for_ER, so not
#      the exact test's rows); unset, or set below every p, nothing is written;
#   3. each dump replays: the kernel named on its second line, fed the vector
#      and the carrier flag, returns the golden row's p; NAmu, NAsigma and var2
#      on that line equal their definitions from the rows (one quantity two
#      ways); the carrier flag count is the golden row's carrier count.
# The first replay found LEDGER #173 (2026-09-28): on the marker path under the
# default --is_noadjCov=TRUE the score test runs on the CENTRED genotype and
# the saddlepoint's CGF on the covariate-PROJECTED one, so the fast kernel's
# Gaussian remainder NAsigma = var2(centred) - carriers(projected) is not the
# non-carriers' variance (P1: 0.91 vs 0.78). Characterized here; the region
# path builds its context with noadjCov = false and is not affected.
skip_if_no_goldens()

marker_args <- function(out) {
  c("--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first",
    paste0("--SAIGEOutputFile=", out),
    "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    paste0("--GMMATmodelFile=", file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
    paste0("--varianceRatioFile=", file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt")),
    "--is_output_moreDetails=TRUE")
}

golden_p1 <- file.path(saige_goldens, "step2", "plink_marker_binary.txt")

# rows whose saddlepoint was attempted and converged: Is.SPA true and MAC above
# the default --max_MAC_for_ER of 4 (LEDGER #80: the exact test's rows read
# Is.SPA true too; they never enter the saddlepoint)
spa_rows_of <- function(gold) {
  ac <- as.numeric(gold$AC_Allele2)
  n2 <- 2 * (as.numeric(gold$N_case) + as.numeric(gold$N_ctrl))
  mac <- pmin(ac, n2 - ac)
  which(tolower(gold$Is.SPA) == "true" & mac > 4)
}

# the dump reader and the replay live in helper-spa-dump.R
read_dump <- read_spa_dump
num <- spa_dump_num
replay <- replay_spa_dump

run_p1_with <- function(env) {
  dir <- tempfile("spadump_"); dir.create(dir)
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", marker_args(out), env = c(SAIGE_SPA_DUMP_DIR = dir, env))
  expect_cli_ok(st)
  list(dir = dir, out = out, log = readLines(attr(st, "log"), warn = FALSE))
}

title_replay <- paste("SAIGE_SPA_DUMP_BELOW=1 dumps every converged saddlepoint of P1,",
                      "moves no byte, and each dump replays to the golden p")
test_that(title_replay, {
  r <- run_p1_with(c(SAIGE_SPA_DUMP_BELOW = "1"))
  expect_output_matches_golden(golden_p1, r$out)
  fallback_line <- grep("^SPA fallbacks: ", r$log, value = TRUE)
  expect_identical(as.integer(sub("^SPA fallbacks: ([0-9]+) .*$", "\\1", fallback_line)), 0L)
  gold <- read_saige_table(golden_p1)
  spa_rows <- spa_rows_of(gold)
  expect_gt(length(spa_rows), 0)
  files <- list.files(r$dir, pattern = "^spaBelow_[0-9]+_[0-9]+\\.tsv$", full.names = TRUE)
  expect_length(list.files(r$dir, pattern = "^spaFail_"), 0L)
  expect_length(files, length(spa_rows))
  matched <- integer(0)
  for (f in files) {
    dm <- read_dump(f)
    expect_identical(dm$h1[["why"]], "converged")
    expect_true(dm$h2[["kernel"]] %in% c("fast", "dense"))
    expect_identical(names(dm$d), c("gtilde", "mu", "y", "nz", "g", "dropped"))
    expect_true(all(dm$d$dropped == 0))   # LEDGER #174: P1's model has no pinned sample
    n_model <- as.integer(as.numeric(gold$N_case[1]) + as.numeric(gold$N_ctrl[1]))
    expect_identical(nrow(dm$d), n_model)
    # the golden row this dump belongs to: Tstat AND var, at the golden's printed precision
    ts <- num(dm$h1, "Tstat")
    v1 <- num(dm$h1, "var1")
    same_tstat <- abs(as.numeric(gold$Tstat) - ts) <= 1e-5 * abs(ts)
    same_var <- abs(as.numeric(gold$var) - v1) <= 1e-5 * abs(v1)
    i <- which(same_tstat & same_var)
    expect_length(i, 1L)
    expect_true(i %in% spa_rows)
    matched <- c(matched, i)
    # the reported p on line 2 is the golden's, and the replay returns it
    p_gold <- as.numeric(gold$p.value[i])
    expect_equal(num(dm$h2, "p"), p_gold, tolerance = 1e-5)
    expect_equal(num(dm$h2, "pNA"), as.numeric(gold$p.value.NA[i]), tolerance = 1e-5)
    rp <- replay(dm)
    expect_true(isTRUE(rp$Isconverge))
    expect_equal(rp$pvalue, p_gold, tolerance = 1e-5)
    # one quantity two ways: the Gaussian remainder from the rows vs the line
    mu <- dm$d$mu
    g <- dm$d$gtilde
    nz <- dm$d$nz == 1
    counts <- c(gold$N_case_hom[i], gold$N_case_het[i], gold$N_ctrl_hom[i], gold$N_ctrl_het[i])
    carriers_gold <- sum(as.numeric(counts))
    expect_equal(sum(nz), carriers_gold)
    fast <- identical(dm$h2[["kernel"]], "fast")
    expect_identical(as.integer(num(dm$h2, "nNB")), if (fast) sum(nz) else 0L)
    m1 <- sum(mu * g)
    expect_equal(num(dm$h1, "m1"), m1, tolerance = 1e-8)
    # the marker path under the default --is_noadjCov=TRUE scores the CENTRED
    # genotype (scoreTestFast_noadjCov): var2 = sum(mu2 * (G - 2f)^2) with f the
    # frequency of the counted allele, NOT the projected vector's variance
    # (LEDGER #173, characterized below)
    G <- dm$d$g
    expect_equal(sum(G != 0), sum(nz))
    f <- sum(G) / (2 * length(G))
    var2_centred <- sum(mu * (1 - mu) * (G - 2 * f)^2)
    expect_equal(num(dm$h2, "var2"), var2_centred, tolerance = 1e-8)
    expect_equal(num(dm$h1, "var1"), num(dm$h2, "var2") * num(dm$h2, "varRatio"), tolerance = 1e-8)
    if (fast) {
      expect_equal(num(dm$h2, "NAmu"), m1 - sum(g[nz] * mu[nz]), tolerance = 1e-8)
      carriers_var <- sum(mu[nz] * (1 - mu[nz]) * g[nz]^2)
      expect_equal(num(dm$h2, "NAsigma"), num(dm$h2, "var2") - carriers_var, tolerance = 1e-6)
    }
    q_expected <- num(dm$h1, "Tstat") / sqrt(num(dm$h2, "varRatio")) + m1
    expect_equal(num(dm$h1, "q"), q_expected, tolerance = 1e-8)
    # LEDGER #173 characterization: the Gaussian remainder the fast kernel was
    # given is NOT the projected vector's non-carrier variance -- the score's
    # var2 is the centred genotype's. Invert when #173 is fixed.
    if (fast) {
      remainder_projected <- sum(mu[!nz] * (1 - mu[!nz]) * g[!nz]^2)
      expect_false(isTRUE(all.equal(num(dm$h2, "NAsigma"), remainder_projected, tolerance = 1e-3)))
      var2_projected <- sum(mu * (1 - mu) * g^2)
      expect_false(isTRUE(all.equal(num(dm$h2, "var2"), var2_projected, tolerance = 1e-3)))
    }
  }
  expect_setequal(matched, spa_rows)
})

test_that("with the threshold unset, or below every p, no converged saddlepoint is dumped", {
  r0 <- run_p1_with(character())
  expect_output_matches_golden(golden_p1, r0$out)
  expect_length(list.files(r0$dir), 0L)
  r1 <- run_p1_with(c(SAIGE_SPA_DUMP_BELOW = "1e-300"))
  expect_output_matches_golden(golden_p1, r1$out)
  expect_length(list.files(r1$dir), 0L)
})
