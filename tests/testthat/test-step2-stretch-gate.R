# LEDGER #176 (2026-10-02): the stretch gate and its shadow kernel.
# OBSERVATION ONLY -- the build this file gates changes no p-value.
#
# The binary saddlepoint is evaluated at q = Tstat / sqrt(varianceRatio) + m1
# against a CGF on the unscaled vector, so the variance ratio stretches the
# observed score. The gate is that stretch over the sd of everything but the
# heaviest carrier's term; the shadow kernel (src/ExactConv.hpp) is the exact
# two-sided tail at the UNSTRETCHED score, by convolution.
#
#   1. The kernel against closed forms: a binomial (mid-p, no Gaussian part),
#      2^n enumeration with a Gaussian part, a null control (q at the mean),
#      and a positive control (a one-atom vector, where a 1 % stretch costs
#      orders of magnitude and a well-spread vector costs almost nothing).
#   2. P78/P77 (the above-cap fixture): the goldens do not move, the run
#      REPORTS the gate, and every value was predicted from the run's own
#      saddlepoint dumps BEFORE the C++ existed (harness/tools/spa_exact_ref.py
#      over SAIGE_SPA_DUMP_BELOW=1 dumps of the unpatched build).
#   3. The floor and the kill switch.

conv <- function(g, mu, rm, rs, q, qinv, grid = 200000L) {
  SAIGE:::exactConvTwoSidedInCPP(g, mu, rm, rs, q, qinv, as.integer(grid))
}

test_that("the convolution reproduces a binomial's mid-p tail (no Gaussian part)", {
  n <- 40L
  p <- 0.1122
  for (k in c(8L, 12L, 20L)) {
    m <- n * p
    qinv <- 2 * m - k
    upper <- pbinom(k, n, p, lower.tail = FALSE) + 0.5 * dbinom(k, n, p)
    lower <- if (qinv < 0) 0 else pbinom(floor(qinv), n, p)   # qinv is off the lattice here
    got <- conv(rep(1, n), rep(p, n), 0, 0, k, qinv, grid = 40L * 100L)
    expect_equal(got$pval, upper + lower, tolerance = 1e-10)
  }
})

test_that("the convolution agrees with 2^n enumeration when there is a Gaussian part", {
  set.seed(11)
  for (trial in 1:6) {
    n <- 10L
    g <- if (trial %% 2 == 0) c(1, runif(n - 1) * 0.05) else rnorm(n)
    mu <- 10^runif(n, -3, -0.4)
    rs <- c(0.05, 0.002)[1 + trial %% 2]
    rm <- rnorm(1, 0, 0.01)
    m1 <- sum(g * mu) + rm
    q <- sum(g[c(1, 3)]) + rm
    if (q < m1) q <- 2 * m1 - q
    qinv <- 2 * m1 - q
    cfg <- as.matrix(expand.grid(rep(list(0:1), n)))
    pr <- apply(cfg, 1, function(b) prod(ifelse(b == 1, mu, 1 - mu)))
    s <- as.vector(cfg %*% g) + rm
    exact <- sum(pr * (pnorm((q - s) / rs, lower.tail = FALSE) + pnorm((qinv - s) / rs)))
    got <- conv(g, mu, rm, rs, q, qinv)
    expect_equal(got$pval, exact, tolerance = 5e-3)
    expect_lte(got$lo, exact * (1 + 1e-9))
    expect_gte(got$hi, exact * (1 - 1e-9))
  }
})

test_that("null control: at the mean the two tails sum to 1", {
  set.seed(3)
  g <- runif(60)
  mu <- runif(60, 0.02, 0.4)
  m1 <- sum(g * mu)
  expect_equal(conv(g, mu, 0, 0.01, m1, m1)$pval, 1, tolerance = 1e-6)
})

test_that("positive control: a 1 % stretch costs orders on one atom and little when spread", {
  vr <- 0.98
  tail_ratio <- function(secondary) {
    g <- c(1, rep(secondary / 81, 81))
    mu <- c(0.01, rep(0.05, 81))
    rm <- -sum(g * mu)                       # a centred score
    obs <- 1 - sum(g * mu)                   # the atom fired, everything else at its mean
    at <- function(t) conv(g, mu, rm, 0.001, t, -t)$pval
    at(obs) / at(obs / sqrt(vr))
  }
  expect_gt(tail_ratio(0.003), 1e6)          # nothing to absorb the stretch
  expect_lt(tail_ratio(2.3), 1.5)            # absorbed
})

# ---- end to end, on the above-cap fixture (goldens P77 and P78) -------------
skip_if_no_goldens()

gate_common <- c("--bedFile=./input/genotype_100markers.bed",
                 "--bimFile=./input/genotype_100markers.bim",
                 "--famFile=./input/genotype_100markers.fam",
                 "--AlleleOrder=alt-first",
                 "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                 "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
                 "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
                 paste0("--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
                        "randomMarkersUsed.sparseGRM.mtx"),
                 paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
                        "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
                 "--is_fastTest=FALSE",
                 "--groupFile=./input/group_new_snpid_withWeights_abovecap.txt",
                 "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
                 "--maxMAF_in_groupTest=0.0001,0.001,0.01")
gate_args <- function(out, burden) {
  c(gate_common, paste0("--SAIGEOutputFile=", out), if (burden) "--r.corr=1" else "--r.corr=0")
}
gate_gold <- function(burden) {
  file.path(saige_goldens, "step2",
            if (burden) "plink_group_customweight_abovecap_burden.txt"
            else        "plink_group_customweight_abovecap_skato.txt")
}
gate_line <- function(log) grep("^Stretch gate: ", log, value = TRUE)
gate_bands <- function(log) {
  as.numeric(regmatches(gate_line(log), gregexpr("[0-9.]+(?= (weighted|at or above))",
                                                 gate_line(log), perl = TRUE))[[1]])
}
gate_side <- function(out) {
  utils::read.table(paste0(out, ".stretchGate.txt"), header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)
}

test_that("P78 in shadow mode: the golden does not move and every value is the predicted one", {
  skip_unless_level("standard")
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", gate_args(out, TRUE),
                      env = c(SAIGE_EXACT_CONV_SHADOW = "1"))
  expect_cli_ok(st)
  log <- readLines(attr(st, "log"), warn = FALSE)
  expect_output_matches_golden(gate_gold(TRUE), out)

  # the run REPORTS it, and the line agrees with the sidecar it wrote
  expect_length(gate_line(log), 1L)
  side <- gate_side(out)
  bands <- gate_bands(log)
  expect_identical(bands[1], as.numeric(nrow(side)))      # shadow mode lists every evaluated cell
  expect_identical(bands[2], as.numeric(sum(side$gate >= 0.5)))
  expect_true(all(side$Kind == "burden"))

  # predicted from the unpatched build's dumps: (Region, max_MAF, first Group)
  want <- data.frame(
    Region  = c("GENE1", "GENE1", "GENE2", "GENE2"),
    Group   = c("lof", "missense;lof", "missense;lof", "missense;lof"),
    max_MAF = c(0.01, 0.01, 0.001, 0.01),
    gate    = c(0.840293, 0.76505, 2.3486e-05, 0.0901574),
    p_exact = c(0.0540288, 0.0527064, 0.0194173, 0.0160107)
  )
  res <- read_saige_table(gate_gold(TRUE))
  res <- res[res$Group != "Cauchy", ]                     # its max_MAF is NA
  for (i in seq_len(nrow(want))) {
    row <- side[side$Region == want$Region[i] & side$Group == want$Group[i] &
                  abs(side$max_MAF - want$max_MAF[i]) < 1e-12, ]
    expect_identical(nrow(row), 1L)
    if (want$gate[i] < 1e-3) {
      expect_lt(abs(row$gate - want$gate[i]), 1e-3)          # a gate of ~0: absolute
    } else {
      expect_equal(row$gate, want$gate[i], tolerance = 1e-3)
    }
    expect_equal(row$p_exact, want$p_exact[i], tolerance = 1e-3)
    expect_lte(row$p_exact_lo, row$p_exact * (1 + 1e-6))
    expect_gte(row$p_exact_hi, row$p_exact * (1 - 1e-6))
    # one statistic two ways: the saddlepoint p the gate recorded is the p the
    # results file printed for that cell (6 printed digits)
    r <- res[res$Region == want$Region[i] & res$Group == want$Group[i] &
               abs(as.numeric(res$max_MAF) - want$max_MAF[i]) < 1e-12, ]
    expect_equal(row$p_spa, as.numeric(r$Pvalue_Burden), tolerance = 1e-5)
  }
  # the fixture has a cell on EACH side of the 0.5 line, and the two above it
  # are the ones whose exact p is far from the saddlepoint's
  over <- side$gate >= 0.5
  expect_true(any(over) && any(!over))
  expect_true(all(side$p_exact[over] / side$p_spa[over] > 1.5))
  expect_true(all(side$p_exact[!over] / side$p_spa[!over] < 1.5))
})

test_that("without shadow mode only the cells at or above the floor are listed", {
  skip_unless_level("standard")
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", gate_args(out, TRUE))
  expect_cli_ok(st)
  log <- readLines(attr(st, "log"), warn = FALSE)
  expect_output_matches_golden(gate_gold(TRUE), out)
  side <- gate_side(out)
  bands <- gate_bands(log)
  expect_true(all(side$gate >= 0.1))
  expect_identical(as.numeric(nrow(side)), bands[3])       # the ">= 0.1" band
  expect_false("p_exact" %in% names(side))
  expect_gt(bands[1], nrow(side))        # some evaluated cells are below the floor
})

test_that("P77: the ultra-rare pass is gated too, and its golden does not move", {
  skip_unless_level("standard")
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", gate_args(out, FALSE),
                      env = c(SAIGE_EXACT_CONV_SHADOW = "1"))
  expect_cli_ok(st)
  expect_output_matches_golden(gate_gold(FALSE), out)
  side <- gate_side(out)
  expect_true(all(side$Kind == "ultra_rare"))
  row <- side[side$Region == "GENE1" & side$Group == "lof" & abs(side$max_MAF - 0.01) < 1e-12, ]
  expect_identical(nrow(row), 1L)
  expect_equal(row$gate, 1.07856, tolerance = 1e-3)        # predicted from the dumps
  expect_equal(row$p_exact, 0.046482, tolerance = 1e-3)
})

test_that("the kill switch silences the gate and writes nothing", {
  skip_unless_level("standard")
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", gate_args(out, TRUE),
                      env = c(SAIGE_STRETCH_GATE = "0"))
  expect_cli_ok(st)
  log <- readLines(attr(st, "log"), warn = FALSE)
  expect_output_matches_golden(gate_gold(TRUE), out)
  expect_match(gate_line(log), "OFF")
  expect_false(file.exists(paste0(out, ".stretchGate.txt")))
})
