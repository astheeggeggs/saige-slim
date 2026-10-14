# --is_output_log10P (2026-10-14): log10 p-value columns on both region paths.
#
# A p-value that underflows double precision (below ~1e-308) printed as 0 is
# a statement about the representation, not the evidence, and a pooled row
# over such a cell was 0 by CCT's `is.zero` rule. The flag adds log10 twins
# of the p-value columns after SE_Burden, exact in log space wherever the
# kernel works there (the SKAT-O saddlepoint tier, the C++ score test's
# underflow spelling, the Cauchy combination), log10(p) elsewhere.
#
# No bundled cell underflows, so the end-to-end runs here pin the arithmetic
# where log10(p) is representable and the goldens P79/P80 pin the bytes; the
# underflow branches are unit-tested on synthetic inputs:
#   (1) off: a run with the flag is the run without it in every pre-existing
#       column, both paths (the flag is additive), and the run REPORTS the flag;
#   (2) on: the columns sit where the schema says and agree with log10 of the
#       printed p on every row including the pooled row (1e-12 on the SKAT-O
#       path, which prints 15 digits; the BURDEN path prints %.6E, so 5e-7);
#   (3) CCT_log10 (R) and CCT_log10_cpp agree with a cancellation-free
#       reference on random vectors INCLUDING underflowed entries, and with
#       the closed form k * p_min when one entry dominates; CCT() itself is
#       measured off that reference (its 1 - pcauchy() cancels at large
#       statistics), which is why the log version does not compare to it;
#   (4) the saddlepoint tier's log twin is finite where the linear p is 0 and
#       equals log(p) where it is not;
#   (5) the C++ digit parser reads the score test's underflow spelling.

log10p_run <- function(path, flag, out = tempfile(fileext = ".txt")) {
  m <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
  sg_name <- "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  sg <- file.path(saige_extdata, "output", sg_name)
  common <- c("--bedFile=./input/genotype_100markers.bed",
              "--bimFile=./input/genotype_100markers.bim",
              "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
              "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
              paste0("--GMMATmodelFile=", m, ".rda"),
              paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
              paste0("--sparseGRMFile=", sg),
              paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
              paste0("--is_output_log10P=", flag),
              paste0("--SAIGEOutputFile=", out))
  args <- if (path == "skato") {
    c(common, "--groupFile=./input/group_new_snpid.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=TRUE")
  } else {
    c(common, "--groupFile=./input/group_new_snpid_flexrv.txt",
      "--flexRV_score=REVEL", "--flexRV_maxMAF=0.01",
      "--annotation_in_groupTest=lof:missense", "--r.corr=1", "--is_fastTest=FALSE")
  }
  st <- run_saige_cli("step2_SPAtests.R", args)
  expect_cli_ok(st)
  list(out = out, log = readLines(attr(st, "log"), warn = FALSE))
}

read_rows <- function(f) read.delim(f, check.names = FALSE, colClasses = "character")

test_that("the flag is additive: every pre-existing column is unchanged, and the run reports it", {
  skip_unless_level("standard")
  for (path in c("skato", "burden")) {
    off <- log10p_run(path, "FALSE")
    on <- log10p_run(path, "TRUE")
    a <- read_rows(on$out)
    b <- read_rows(off$out)
    new_cols <- grep("^log10_", names(a), value = TRUE)
    expect_gt(length(new_cols), 0)
    expect_identical(a[, setdiff(names(a), new_cols)], b, info = path)
    expect_false(any(grepl("^log10_", names(b))), info = path)
    expect_true(any(grepl("log10 p-value columns (--is_output_log10P) are on", on$log,
                          fixed = TRUE)), info = path)
    expect_false(any(grepl("--is_output_log10P) are on", off$log, fixed = TRUE)), info = path)
  }
})

test_that("SKAT-O path: three columns after SE_Burden equal log10(printed p), pooled row too", {
  skip_unless_level("standard")
  d <- read_rows(log10p_run("skato", "TRUE")$out)
  i <- match("SE_Burden", names(d))
  expect_identical(names(d)[i + 1:3], c("log10_Pvalue", "log10_Pvalue_Burden", "log10_Pvalue_SKAT"))
  for (pr in list(c("Pvalue", "log10_Pvalue"), c("Pvalue_Burden", "log10_Pvalue_Burden"),
                  c("Pvalue_SKAT", "log10_Pvalue_SKAT"))) {
    p <- suppressWarnings(as.numeric(d[[pr[1]]]))
    l <- suppressWarnings(as.numeric(d[[pr[2]]]))
    expect_identical(is.na(p), is.na(l), info = pr[2])
    ok <- !is.na(p)
    expect_gt(sum(ok), 0)
    expect_lt(max(abs(log10(p[ok]) - l[ok])), 1e-12)
  }
  expect_true(any(d$Group == "Cauchy" & !is.na(suppressWarnings(as.numeric(d$log10_Pvalue)))))
})

test_that("BURDEN path: the column follows SE_Burden and equals log10(printed p), pooled row too", {
  skip_unless_level("standard")
  d <- read_rows(log10p_run("burden", "TRUE")$out)
  i <- match("SE_Burden", names(d))
  expect_identical(names(d)[i + 1], "log10_Pvalue_Burden")
  p <- suppressWarnings(as.numeric(d$Pvalue_Burden))
  l <- suppressWarnings(as.numeric(d$log10_Pvalue_Burden))
  expect_identical(is.na(p), is.na(l))
  ok <- !is.na(p)
  expect_gt(sum(ok), 100)
  # a cell's p is printed %.6E and its log10 column is log10 of those printed
  # digits, so the two agree to the column's own %.10f rounding; the pooled
  # row's p is printed by writeNumOrNA (fewer digits) while its log10 column is
  # the log-space combination of the cells' digits, so that row agrees only to
  # the printed p's precision
  cell <- ok & d$Group != "Cauchy"
  expect_lt(max(abs(log10(p[cell]) - l[cell])), 1e-9)
  pooled <- ok & d$Group == "Cauchy"
  expect_gt(sum(pooled), 0)
  expect_lt(max(abs(log10(p[pooled]) - l[pooled])), 5e-6)
  expect_true(any(d$Group == "Cauchy" & !is.na(l)))
})

# a cancellation-free reference for the Cauchy combination's log10 p
cct_log10_ref <- function(p) {
  k <- length(p)
  sm <- p < 1e-16
  st <- (sum(1 / (p[sm] * pi)) + sum(tan((0.5 - p[!sm]) * pi))) / k
  if (st > 1e15) return(-(log10(st) + log10(pi)))
  pcauchy(st, lower.tail = FALSE, log.p = TRUE) / log(10)
}

test_that("CCT_log10 and CCT_log10_cpp match a cancellation-free reference, underflow included", {
  set.seed(20261014)
  worst_r <- 0
  worst_c <- 0
  worst_cct <- 0
  for (i in seq_len(1500)) {
    k <- sample(1:30, 1)
    p <- 10^runif(k, -14, 0)
    if (runif(1) < 0.3) p[sample(k, 1)] <- 10^runif(1, -300, -17)
    r <- cct_log10_ref(p)
    a <- SAIGE:::CCT_log10(log10(p))
    cc <- CCT_log10_cpp(log10(p))
    worst_r <- max(worst_r, abs(a - r))
    worst_c <- max(worst_c, abs(cc - a))
    worst_cct <- max(worst_cct, abs(log10(SAIGE:::CCT(p)) - r))
  }
  expect_lt(worst_r, 1e-10)
  expect_lt(worst_c, 1e-12)
  # the control that the reference is needed: CCT()'s own 1 - pcauchy() loses
  # digits at large statistics (measured 2e-3 in log10 p on this draw)
  expect_gt(worst_cct, 1e-6)
  # one entry dominates: the combination is k * p_min, exactly, even where p_min
  # is far below double precision
  x <- c(-400, -2, -1.3, -0.5)
  expect_equal(SAIGE:::CCT_log10(x), -400 + log10(4), tolerance = 1e-12)
  expect_equal(CCT_log10_cpp(x), -400 + log10(4), tolerance = 1e-12)
  # a p of 1 anywhere: CCT()'s min(1, k * min p) rule
  expect_equal(SAIGE:::CCT_log10(c(-3, 0)), log10(2e-3), tolerance = 1e-12)
  expect_equal(CCT_log10_cpp(c(-3, 0)), log10(2e-3), tolerance = 1e-12)
  # no inputs
  expect_true(is.na(SAIGE:::get_CCT_log10pvalue(c(NA_real_, NA_real_))))
})

test_that("the saddlepoint log twin is finite where the linear p underflows, log(p) elsewhere", {
  lambda <- c(2.5, 1.2, 0.7, 0.3, 0.1)
  q_mod <- 200                                  # representable: p ~ 7.6e-19
  p <- SAIGE:::skatoSaddlepointSurvival(lambda, q_mod)
  lp <- SAIGE:::skatoSaddlepointSurvival(lambda, q_mod, log = TRUE)
  expect_true(p > 0 && p < 1e-15)
  expect_equal(lp, log(p), tolerance = 1e-12)
  q_big <- 5000                                 # the linear value underflows
  p2 <- SAIGE:::skatoSaddlepointSurvival(lambda, q_big)
  lp2 <- SAIGE:::skatoSaddlepointSurvival(lambda, q_big, log = TRUE)
  expect_identical(p2, 0)
  expect_true(is.finite(lp2) && lp2 < -700)
  expect_identical(SAIGE:::skatoSaddlepointSurvival(lambda, 0, log = TRUE), 0)
})

test_that("the C++ digit parser reads the score test's underflow spelling", {
  expect_equal(log10PrintedPval_cpp("2.500000E-03"), log10(2.5e-3), tolerance = 1e-12)
  # formatLogP's spelling, past DBL_MIN
  expect_equal(log10PrintedPval_cpp("1.0E-400"), -400, tolerance = 1e-12)
  expect_equal(log10PrintedPval_cpp("3.7E-320"), log10(3.7) - 320, tolerance = 1e-12)
  expect_equal(log10PrintedPval_cpp("0.5"), log10(0.5), tolerance = 1e-12)
  expect_true(is.na(log10PrintedPval_cpp("NA")))
  expect_true(is.na(log10PrintedPval_cpp("0.000000E+00")))
})
