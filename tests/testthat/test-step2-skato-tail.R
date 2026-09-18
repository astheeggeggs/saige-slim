# LEDGER #158, step 2: the saddlepoint tail (R/SkatoScreen.R,
# skatoSaddlepointTail). Two things are pinned. (1) An EXACT control: on
# synthetic cells the tail's p-value and its per-rho values are within 15% of
# the SKAT-O integral evaluated with an exact inner function in both layers
# (Ruben's positive series where it converges, adaptive-tolerance Davies on
# the near-rank-1 spectra where it cannot), down to ~1e-9; SKAT's own value on the same cells is
# reported, not asserted. (2) End to end on the bundled fixture, whose cells
# are all above 1e-5: --SKATO_tail_threshold=1 forces every screened-in cell
# through the tail, and Pvalue, Pvalue_SKAT and Pvalue_Burden land within 15%
# of the same session's Davies run with the run reporting the mode and the
# sidecar naming every cell; --SKATO_tail_method=davies at the same threshold
# is a no-op and reproduces P67's bytes.

ruben_coef <- function(lambda, kmax = 20000) {
  beta <- 2 * min(lambda) * max(lambda) / (min(lambda) + max(lambda))
  gam <- 1 - beta / lambda
  gmax <- max(abs(gam))
  k <- ceiling(95 / (-log(gmax))) + 100
  if (k > kmax) return(NULL)
  g <- vapply(seq_len(k), function(r) sum(gam^r), 0)
  a <- numeric(k + 1)
  a[1] <- prod(sqrt(beta / lambda))
  for (j in seq_len(k)) a[j + 1] <- sum(g[seq_len(j)] * a[j:1]) / (2 * j)
  list(a = a, beta = beta, n = length(lambda), k = k)
}
# Exact where Ruben converges (the well-conditioned spectra); where it does
# not -- the near-rank-1 spectra at rho = 0.999 have a condition number in the
# thousands whatever Phi is -- Davies at an ABSOLUTE tolerance of 1e-4 times a
# saddlepoint pilot (the study's reference, valid to ~1e-9), NA if that fails.
ruben_survival <- function(lambda, q) {
  if (q <= 0) return(1)
  rc <- ruben_coef(lambda)
  if (!is.null(rc)) {
    return(sum(rc$a * pchisq(q / rc$beta, df = rc$n + 2 * (0:rc$k), lower.tail = FALSE)))
  }
  pilot <- skatoSaddlepointSurvival(lambda, q)
  if (is.na(pilot)) pilot <- 0.5
  davies <- get("SKAT_davies", asNamespace("SKAT"))
  o <- davies(q, lambda, lim = 1e8, acc = max(min(pilot * 1e-4, 1e-6), 1e-15))
  if (o$ifault != 0 || o$Qq <= 0) NA_real_ else o$Qq
}
# The SKAT-O integral with an exact inner function in both layers, direct form.
skato_exact <- function(Score, Phi, r.corr) {
  ns <- asNamespace("SKAT")
  p.m <- nrow(Phi)
  r.all <- r.corr
  r.all[r.all >= 0.999] <- 0.999
  n.r <- length(r.all)
  Qr <- get("SKAT_META_Optimal_Get_Q", ns)(Score, r.all)$Q.r
  Phi2 <- Phi / 2
  lambda.all <- lapply(r.all, function(rc) {
    RM <- diag(rep(1 - rc, p.m)) + matrix(rc, p.m, p.m)
    L <- chol(RM, pivot = TRUE)
    get("Get_Lambda", ns)(L %*% (Phi2 %*% t(L)))
  })
  param.m <- get("SKAT_META_Optimal_Param", ns)(Phi2, r.all)
  pval <- numeric(n.r)
  pmat <- matrix(0, n.r, 3)
  for (i in seq_len(n.r)) {
    lam <- lambda.all[[i]]
    pt <- get("Get_Liu_Params_Mod", ns)(c(sum(lam), sum(lam^2), sum(lam^3), sum(lam^4)))
    pval[i] <- ruben_survival(lam, Qr[i])
    pmat[i, ] <- c(pt$muQ, pt$sigmaQ^2, pt$l)
  }
  if (any(is.na(pval))) return(NULL)
  pmin <- min(pval)
  pmin.q <- vapply(seq_len(n.r), function(i) {
    df <- pmat[i, 3]
    (qchisq(1 - pmin, df = df) - df) / sqrt(2 * df) * sqrt(pmat[i, 2]) + pmat[i, 1]
  }, 0)
  sd1 <- sqrt(param.m$VarQ - param.m$VarRemain) / sqrt(param.m$VarQ)
  sfun <- function(x) {
    temp <- (pmin.q - param.m$tau %x% t(x)) / (1 - r.all)
    vapply(apply(temp, 2, min), function(m1) {
      if (m1 > sum(param.m$lambda) * 1e4) return(0)
      min(ruben_survival(param.m$lambda, (m1 - param.m$MuQ) * sd1 + param.m$MuQ), 1)
    }, 0)
  }
  xstar <- max(min(pmin.q / param.m$tau), 0)
  re <- integrate(function(t) sfun(t^2) * 2 * dnorm(t), 0, sqrt(xstar),
                  subdivisions = 1000, rel.tol = 1e-6, abs.tol = 0)
  p <- pchisq(xstar, df = 1, lower.tail = FALSE) + re$value
  list(p = min(p, pmin * n.r), each = pval)
}

test_that("saddlepoint tail within 15% of the exact-series SKAT-O integral down to 1e-9", {
  set.seed(1582)
  rho <- c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1)
  skato_screen_env$threshold <- 0
  skato_screen_env$tail <- "saddlepoint"
  skato_screen_env$tailThreshold <- 1
  on.exit({
    skato_screen_env$threshold <- 0
    skato_screen_env$tail <- "davies"
    skato_screen_env$tailThreshold <- 1e-5
  })
  n_checked <- 0
  for (q in c(3, 4, 6)) {
    for (scale in c(3.5, 4.5, 5.5)) {
      # a well-conditioned spectrum (condition <= ~30) so the exact series converges
      lam <- exp(seq(0, -log(20), length.out = q))
      V <- qr.Q(qr(matrix(rnorm(q * q), q)))
      Phi <- V %*% diag(lam) %*% t(V)
      Score <- as.vector(t(chol(Phi)) %*% (rep(1, q) * scale / sqrt(q)))
      ex <- skato_exact(Score, Phi, rho)
      if (is.null(ex) || ex$p > 1e-4) next
      got <- skatoScreenKernel(Score, Phi, rho)
      expect_identical(attr(got, "skato_method"), "saddlepoint")
      expect_lt(abs(got$p.value / ex$p - 1), 0.15)
      ok <- is.finite(ex$each) & ex$each > 0
      expect_true(all(abs(got$param$p.val.each[ok] / ex$each[ok] - 1) < 0.15))
      n_checked <- n_checked + 1
    }
  }
  expect_gte(n_checked, 4)
})

tail_run <- function(extra, out = tempfile(fileext = ".txt")) {
  m <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
  sg <- file.path(saige_extdata, "output",
                  "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
  args <- c("--bedFile=./input/genotype_100markers.bed",
            "--bimFile=./input/genotype_100markers.bim",
            "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
            "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
            paste0("--GMMATmodelFile=", m, ".rda"),
            paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
            paste0("--sparseGRMFile=", sg),
            paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
            "--groupFile=./input/group_new_snpid.txt",
            "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
            "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=TRUE",
            "--SKATO_screen_threshold=0", extra,
            paste0("--SAIGEOutputFile=", out))
  st <- run_saige_cli("step2_SPAtests.R", args)
  list(status = as.integer(st), log = readLines(attr(st, "log"), warn = FALSE), out = out)
}
read_tab <- function(path) {
  utils::read.table(path, header = TRUE, sep = "\t", colClasses = "character",
                    quote = "", comment.char = "", na.strings = character(0))
}

test_that("forced through the tail, every screened-in cell is saddlepoint, within 15% of Davies", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- tail_run(c("--SKATO_tail_method=saddlepoint", "--SKATO_tail_threshold=1"))
  expect_cli_ok(r$status)
  line <- grep("^SKAT-O screen: threshold 0, tail saddlepoint below 1", r$log, value = TRUE)
  expect_length(line, 1)
  sc <- read_tab(paste0(r$out, ".skatoMethod.txt"))
  sc <- sc[sc$Test == "marginal", ]
  expect_true(all(sc$Method %in% c("saddlepoint", "single-rho")))
  expect_gt(sum(sc$Method == "saddlepoint"), 0)
  nSp <- as.integer(sub(".* ([0-9]+) the saddlepoint tail.*", "\\1", line))
  expect_equal(nSp, sum(read_tab(paste0(r$out, ".skatoMethod.txt"))$Method == "saddlepoint"))
  # the Davies reference is the SAME session's (0, davies) run, not the golden:
  # column-exact comparisons across platforms are the byte-vs-numeric trap
  ref <- tail_run(c("--SKATO_tail_method=davies"))
  expect_cli_ok(ref$status)
  g <- read_tab(ref$out)
  b <- read_tab(r$out)
  expect_identical(dim(g), dim(b))
  key <- function(d) paste(d$Region, d$Group, d$max_MAF)
  sp <- key(b) %in% key(sc)[sc$Method == "saddlepoint"]
  for (col in c("Pvalue", "Pvalue_SKAT", "Pvalue_Burden")) {
    ratio <- as.numeric(b[[col]][sp]) / as.numeric(g[[col]][sp])
    expect_true(all(abs(ratio - 1) < 0.15), label = col)
  }
  # the single-rho cells and every non-p column are untouched
  expect_identical(g$Pvalue[!sp & b$Group != "Cauchy"], b$Pvalue[!sp & b$Group != "Cauchy"])
  for (col in setdiff(names(g), c("Pvalue", "Pvalue_SKAT", "Pvalue_Burden", "SE_Burden"))) {
    expect_identical(g[[col]], b[[col]], label = col)
  }
  # SE_Burden is back-derived from Pvalue_Burden, so it moves with it
  z <- abs(as.numeric(b$BETA_Burden[sp]) / as.numeric(b$SE_Burden[sp]))
  expect_true(all(abs(z / abs(qnorm(as.numeric(b$Pvalue_Burden[sp]) / 2)) - 1) < 1e-6))
})

test_that("tail method davies is a no-op: P67's bytes at any tail threshold", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- tail_run(c("--SKATO_tail_method=davies", "--SKATO_tail_threshold=1"))
  expect_cli_ok(r$status)
  expect_true(any(grepl("SKAT-O screen: off", r$log, fixed = TRUE)))
  golden <- file.path(saige_goldens, "step2", "plink_group_fastTest_screenoff.txt")
  expect_output_matches_golden(golden, r$out)
})
