# Conditional analysis at NON-NEGLIGIBLE conditioning -- the regime no bundled
# fixture reaches, and the only route to LEDGER #66 and #84.
#
# On every committed fixture the conditioning is nearly a no-op: golden P12's
# conditional SPA clears its gate on 4 of 70 markers with var/var_c between
# 1.000032 and 1.007334, and the best var/var_c among markers clearing even the
# LOWEST reachable gate (--SPAcutoff is floored at 0.5) is 1.008504. So the
# conditional variance subtraction and the conditional SPA have only ever been
# exercised where there is nothing to subtract.
#
# The fixture is DERIVED AT RUN TIME and carries no golden. It is built from a
# seed, its markers sit at designed LD with a conditioning marker, and one is
# correlated with the null model's own residual so that a large conditional
# statistic survives the conditioning. There is deliberately no byte comparison:
# the genotypes come from R's RNG, so a golden would pin the RNG rather than the
# code. Every quantity the assertions use is MEASURED on the analysed samples
# inside the test, so a fixture that failed to build is distinguishable from a
# behaviour change.

ld_seed <- 20260904L
ld_maf <- 0.25

# helper-saige.R supplies these; testthat sources it first. Bound at top level
# so this file is self-contained for static analysis as well as at run time.
ld_src_prefix <- file.path(saige_extdata, "input", "genotype_100markers")
ld_model_prefix <- file.path(saige_goldens, "step1", "binary_sparseGRM_a")
cli_run <- run_saige_cli
cli_ok <- expect_cli_ok
read_tbl <- read_saige_table

# Build bed/bim/fam over the SAME fam rows as genotype_100markers, so the
# committed binary step-1 model can be reused unchanged.
build_ld_fixture <- function(dir) {
  src <- ld_src_prefix
  e <- new.env()
  load(paste0(ld_model_prefix, ".rda"), envir = e)
  res <- as.vector(e$modglmm$residuals)
  fam <- utils::read.table(paste0(src, ".fam"), stringsAsFactors = FALSE)
  n <- nrow(fam)
  inmodel <- match(e$modglmm$sampleID, fam[[2]])
  stopifnot(!anyNA(inmodel))

  set.seed(ld_seed)
  thr <- stats::qnorm(c((1 - ld_maf)^2, 1 - ld_maf^2))
  geno_from <- function(z) as.integer(z > thr[1]) + as.integer(z > thr[2])
  # a fraction f of samples redrawn independently gives cor = 1 - f, r2 = (1-f)^2
  partner <- function(g, r2) {
    i <- which(stats::runif(n) < 1 - sqrt(r2))
    g[i] <- geno_from(stats::rnorm(length(i)))
    g
  }
  zsig <- stats::rnorm(n)
  zsig[inmodel] <- 0.30 * scale(res)[, 1] + sqrt(1 - 0.30^2) * stats::rnorm(length(inmodel))
  zpar <- sqrt(0.5) * zsig + sqrt(0.5) * stats::rnorm(n)

  cond <- geno_from(stats::rnorm(n))
  mk <- list(cond = cond, ld100 = cond,
             ld090 = partner(cond, 0.90), ld050 = partner(cond, 0.50),
             ld000 = geno_from(stats::rnorm(n)),
             sig = geno_from(zsig), sigpart = geno_from(zpar))

  prefix <- file.path(dir, "ldfix")
  writeLines(sprintf("1\t%s\t0\t%d\tA\tT", names(mk), seq_along(mk)),
             paste0(prefix, ".bim"))
  file.copy(paste0(src, ".fam"), paste0(prefix, ".fam"), overwrite = TRUE)
  code <- c(3L, 2L, 0L)
  nb <- ceiling(n / 4)
  con <- file(paste0(prefix, ".bed"), "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  for (g in mk) {
    cd <- c(code[g + 1L], rep(0L, nb * 4L - n))
    m <- matrix(cd, nrow = 4L)
    writeBin(as.raw(m[1, ] + 4L * m[2, ] + 16L * m[3, ] + 64L * m[4, ]), con)
  }
  gs <- lapply(mk, function(g) g[inmodel])
  list(prefix = prefix,
       r2_cond = vapply(gs, function(g) stats::cor(gs$cond, g)^2, numeric(1)),
       r2_sig = vapply(gs, function(g) stats::cor(gs$sig, g)^2, numeric(1)),
       cor_res = vapply(gs, function(g) stats::cor(g, res), numeric(1)))
}

run_cond <- function(prefix, out, noadjcov) {
  m <- ld_model_prefix
  args <- c(paste0("--bedFile=", prefix, ".bed"),
            paste0("--bimFile=", prefix, ".bim"),
            paste0("--famFile=", prefix, ".fam"),
            "--AlleleOrder=alt-first",
            paste0("--SAIGEOutputFile=", out),
            "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
            paste0("--GMMATmodelFile=", m, ".rda"),
            paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
            "--is_output_moreDetails=TRUE", "--condition=sigpart,cond",
            paste0("--is_noadjCov=", ifelse(noadjcov, "TRUE", "FALSE")))
  st <- cli_run("step2_SPAtests.R", args)
  cli_ok(st)
  d <- read_tbl(out)
  rownames(d) <- d$MarkerID
  d
}

num <- function(x) suppressWarnings(as.numeric(x))

test_that("the constructed LD fixture is what it claims to be", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  dir <- tempfile()
  dir.create(dir)
  f <- build_ld_fixture(dir)
  # measured on the 1,000 ANALYSED samples, not on all 10,000 fam rows
  expect_equal(unname(f$r2_cond[["ld100"]]), 1, tolerance = 1e-12)
  expect_gt(f$r2_cond[["ld090"]], 0.85)
  expect_lt(f$r2_cond[["ld090"]], 0.95)
  expect_gt(f$r2_cond[["ld050"]], 0.40)
  expect_lt(f$r2_cond[["ld050"]], 0.60)
  expect_lt(f$r2_cond[["ld000"]], 0.01)              # the negative control's arm
  expect_gt(abs(f$cor_res[["sig"]]), 0.25)           # sig really is associated
  expect_gt(f$r2_sig[["sigpart"]], 0.15)             # and sigpart really is in LD with it
  expect_lt(f$r2_sig[["sigpart"]], 0.60)
})

test_that("LEDGER #84: the conditional SPA runs at non-negligible conditioning", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  dir <- tempfile()
  dir.create(dir)
  f <- build_ld_fixture(dir)
  d <- run_cond(f$prefix, tempfile(fileext = ".txt"), noadjcov = TRUE)

  stat <- function(r) num(d[r, "Tstat"])^2 / num(d[r, "var"])
  stat_c <- function(r) num(d[r, "Tstat_c"])^2 / num(d[r, "var_c"])

  # NEGATIVE CONTROL: conditioning on an independent marker is a near no-op.
  expect_lt(abs(num(d["ld000", "var"]) / num(d["ld000", "var_c"]) - 1), 0.01)
  expect_lt(abs(stat_c("ld000") / stat("ld000") - 1), 0.10)

  # POSITIVE CONTROL: at r2 ~ 0.5 the variance really is reduced, and stays
  # finite and ordinary. 1/(1-r2) is the expected order of magnitude.
  vr50 <- num(d["ld050", "var"]) / num(d["ld050", "var_c"])
  expect_gt(vr50, 1.3)
  expect_lt(vr50, 4)

  # THE POINT: a marker whose conditional statistic clears the SPA gate
  # (SPAcutoff = 2 by default, so stat_c > 4) WHILE the conditioning is real.
  # Nothing in the committed goldens does this -- their best var/var_c among
  # markers clearing even the lowest reachable gate is 1.0085.
  expect_gt(stat_c("sig"), 4)
  vr <- num(d["sig", "var"]) / num(d["sig", "var_c"])
  expect_gt(vr, 1.3)
  # and the conditional SPA actually moved the p-value: the saddlepoint p and
  # the score p differ, which is the direct evidence the block executed.
  expect_false(identical(d["sig", "p.value_c"], d["sig", "p.value.NA_c"]))
  decades <- abs(log10(num(d["sig", "p.value_c"])) -
                   log10(num(d["sig", "p.value.NA_c"])))
  expect_gt(decades, 1)
})

# --- LEDGER #66, FIXED: the conditional variance's relative cancellation guard -
# The conditional variance is a DIFFERENCE, var1 - G1'PG2 VarInv G2'PG1, so a
# marker collinear with the conditioning set cancels it to float noise. The only
# guard used to be `varT_c <= DBL_MIN`, which is ~294 orders of magnitude below
# where that cancellation actually lands -- the same geometry LEDGER #44 had to
# fix in the two "fast" score tests. It now uses #44's relative test, with #44's
# threshold and #44's helper, so a degenerate cell says NA instead of publishing
# a ratio of two cancelled noise terms.
test_that("LEDGER #66 FIXED: collinear conditioning reports NA, not cancelled noise", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  dir <- tempfile()
  dir.create(dir)
  f <- build_ld_fixture(dir)
  adj <- run_cond(f$prefix, tempfile(fileext = ".txt"), noadjcov = FALSE)
  def <- run_cond(f$prefix, tempfile(fileext = ".txt"), noadjcov = TRUE)

  # varT_c = var1 - varT_ctemp, so the guard's term scale is
  # |var1| + |varT_ctemp| = 2*var - var_c, computable from the printed columns
  # alone -- no instrumented build needed.
  ratio <- function(d, r) {
    v <- num(d[r, "var"])
    vc <- num(d[r, "var_c"])
    vc / (2 * v - vc)
  }
  # n*eps at the model's 1,000 samples, which is the threshold the C++ uses
  thr <- 1000 * .Machine$double.eps

  # ld100 is byte-identical to the conditioning marker `cond`, so the true
  # conditional variance is 0.
  for (r in c("cond", "ld100", "sigpart")) {
    vc <- num(adj[r, "var_c"])
    # still cancelled to float noise -- the fix does not change the arithmetic
    expect_lt(vc, 1e-10)
    # ... and still ~294 orders of magnitude ABOVE the absolute guard at
    # src/SAIGE_test.cpp, which is why an absolute test could never catch it
    expect_gt(vc, 1e-100)
    # what the RELATIVE test sees: machine epsilon, far below n*eps
    expect_lt(ratio(adj, r), thr / 10)
    # so the effect size is declined rather than published
    expect_identical(adj[r, "BETA_c"], "NA")
    expect_identical(adj[r, "SE_c"], "NA")
    expect_identical(adj[r, "p.value_c"], "1.000000E+00")
  }

  # POSITIVE CONTROL: the guard is scale-free, not a blanket suppression. The
  # r2 ~ 0.9 partner has a genuinely reduced conditional variance and reports an
  # ordinary finite BETA_c/SE_c -- and its ratio sits ~11 orders above the
  # threshold, so nothing here is tuned to the fixture.
  for (r in c("ld090", "ld050", "ld000", "sig")) {
    expect_gt(ratio(adj, r), thr * 1e8)
    expect_true(is.finite(num(adj[r, "BETA_c"])))
    expect_true(is.finite(num(adj[r, "SE_c"])))
  }

  # CONTRAST, and the reason no CLI-driven golden could have caught this: at the
  # CLI default (--is_noadjCov=TRUE) the tested marker is NOT covariate-adjusted
  # while the conditioning marker's factors always are, so the subtraction does
  # not cancel and the same perfectly collinear marker reports an ORDINARY
  # finite variance instead. The degenerate regime is gated behind the flag
  # whose API default is FALSE and whose CLI default is TRUE (LEDGER #82).
  expect_gt(num(def["ld100", "var_c"]), 1e-3)
  expect_lt(num(def["ld100", "SE_c"]), 1e3)
  expect_gt(ratio(def, "ld100"), thr * 1e8)
})
