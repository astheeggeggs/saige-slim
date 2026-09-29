# LEDGER #174 (2026-09-29): PINNED samples and the saddlepoint.
#
# A sample whose outcome the null model holds certain -- flip probability (mu
# for a control, 1 - mu for a case) at or below 1e-8 -- adds nothing to the
# score or its variance, but its term in the saddlepoint's CGF is inert while
# g t is below -log(mu) and explosive beyond. The binary saddlepoint leaves
# such samples out of its vector while the mass left out of one test is within
# the budget (1e-6). Two models (harness/make_goldens.sh, P72-P76):
#   PLANTED  the separated group's fitted values at 1e-12: the regime the guard
#            was designed on, in a cell LEDGER #172's routing does not take;
#   FITTED   the separated group where step 1 left it, 1.7e-9..1.2e-7: 13 of
#            the 60 are at or below the threshold and 47 are not.
# What is held here:
#   1. the fixture's facts, from the committed models;
#   2. the REFERENCE for each tested cell, two ways (an enumeration of the
#      material carriers against a normal remainder; a seeded Monte Carlo);
#   3. planted: guard off, the saddlepoint is wrong by more than four orders
#      of magnitude on both region paths; guard on (goldens P74, P75), it is
#      within 15 % of the reference; every other row is where it was;
#   4. fitted: guard off, wrong by 1.5x and 2.8x; at the SHIPPED threshold
#      (golden P72) it is still wrong -- a characterization, to be inverted
#      when the threshold is decided; at 1e-6 (golden P73) it is within 15 %;
#   5. the budget: over it nothing is left out and the run says so;
#   6. one quantity two ways: a guarded dump replays to the run's p;
#   7. controls: an unweighted scan does not move (golden P76), a model with
#      no pinned sample reports none, a setting that does not parse is an error.
skip_if_no_goldens()

pin_prefix <- function(planted) {
  file.path(saige_goldens, "step1",
            if (planted) "binary_sparseGRM_pinned_planted" else "binary_sparseGRM_pinned")
}
pin_base <- function(out, planted) {
  c("--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first", "--chrom=1", "--minMAF=0", "--LOCO=FALSE",
    paste0("--SAIGEOutputFile=", out),
    paste0("--GMMATmodelFile=", pin_prefix(planted), ".rda"),
    paste0("--varianceRatioFile=", pin_prefix(planted), ".varianceRatio.txt"))
}
sparse_grm <- c(
  paste0("--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
         "randomMarkersUsed.sparseGRM.mtx"),
  paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
         "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
)
region_args <- function(out, planted, burden) {
  group <- if (planted) "pinnedTail" else "pinned"
  c(pin_base(out, planted), "--minMAC=0.5",
    paste0("--groupFile=./input/group_new_snpid_withWeights_", group, ".txt"),
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=FALSE",
    if (burden) "--r.corr=1" else sparse_grm)
}
marker_args <- function(out, planted) {
  c(pin_base(out, planted), "--minMAC=1", "--is_output_moreDetails=TRUE")
}
gold <- function(f) file.path(saige_goldens, "step2", f)
eps6 <- c(SAIGE_SPA_PINNED_EPS = "1e-6", SAIGE_SPA_PINNED_BUDGET = "1e-5")
guard_off <- c(SAIGE_SPA_PINNED_GUARD = "0")

run_pinned <- function(args_fun, env = character(), dump = FALSE, ...) {
  out <- tempfile(fileext = ".txt")
  dir <- NULL
  if (dump) {
    dir <- tempfile("spadump_")
    dir.create(dir)
    env <- c(env, SAIGE_SPA_DUMP_DIR = dir, SAIGE_SPA_DUMP_BELOW = "1")
  }
  st <- run_saige_cli("step2_SPAtests.R", args_fun(out, ...), env = env)
  expect_cli_ok(st)
  list(out = out, dir = dir, log = readLines(attr(st, "log"), warn = FALSE))
}
pinned_line <- function(log) grep("^Pinned samples: ", log, value = TRUE)
read_side <- function(path) {
  utils::read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}
# the two cells of GENE1 that hold the separated carriers (max MAF 0.01)
cells_01 <- function(tab) which(tab$Region == "GENE1" & tab$max_MAF == "0.01")
pnum <- function(x) as.numeric(x)

# ---- the reference ----------------------------------------------------------
# The samples whose outcome is random at the p-values in play (flip probability
# above 1e-6) are the BODY. The reference is the two-sided tail of the body's
# sum at the run's centred score, computed by enumerating the body's material
# carriers (those within a factor 20 of its largest entry) against a normal
# remainder -- the remainder is hundreds of entries a hundredth of that size.
reference_p <- function(dm) {
  d <- dm$d
  flip <- ifelse(d$y == 1, 1 - d$mu, d$mu)
  body <- flip > 1e-6
  g <- d$gtilde[body]
  mu <- d$mu[body]
  centred <- abs(as.numeric(dm$h1[["q"]]) - as.numeric(dm$h1[["m1"]]))
  big <- abs(g) >= 0.05 * max(abs(g))
  stopifnot(sum(big) >= 1, sum(big) <= 16)
  gm <- g[big]
  mm <- mu[big]
  m_rest <- sum(mu[!big] * g[!big])
  sd_rest <- sqrt(sum(mu[!big] * (1 - mu[!big]) * g[!big]^2))
  m_body <- sum(mu * g)
  p <- 0
  for (cfg in 0:(2^length(gm) - 1)) {
    y <- as.integer(intToBits(cfg))[seq_along(gm)]
    s <- sum(gm * y)
    pr <- prod(ifelse(y == 1, mm, 1 - mm))
    p <- p + pr * (stats::pnorm(m_body + centred - s, m_rest, sd_rest, lower.tail = FALSE) +
                     stats::pnorm(m_body - centred - s, m_rest, sd_rest))
  }
  list(p = p, body = body, centred = centred, m_body = m_body, n_material = sum(big))
}
monte_carlo_p <- function(dm, ref, draws = 1e5) {
  g <- dm$d$gtilde[ref$body]
  mu <- dm$d$mu[ref$body]
  set.seed(174)
  hits <- 0
  chunk <- 1e4
  for (k in seq_len(draws / chunk)) {
    y <- matrix(stats::rbinom(length(g) * chunk, 1, rep(mu, times = chunk)), nrow = length(g))
    s <- as.vector(crossprod(g, y))
    hits <- hits + sum(s >= ref$m_body + ref$centred | s <= ref$m_body - ref$centred)
  }
  c(p = hits / draws, se = sqrt(hits) / draws)
}
# the dumps of one run, keyed by the cell they belong to (largest carrier count
# is the larger cell: missense;lof holds lof's markers and more)
dumps_of <- function(dir) {
  files <- list.files(dir, pattern = "^spaBelow_[0-9]+_[0-9]+\\.tsv$", full.names = TRUE)
  lapply(files, read_spa_dump) # nolint: object_usage_linter.
}
within <- function(p, ref, frac) abs(p / ref - 1) <= frac

# ---- 1. the fixture ---------------------------------------------------------
load_model <- function(planted) {
  e <- new.env()
  load(paste0(pin_prefix(planted), ".rda"), envir = e)
  e$modglmm
}

title_fixture <- paste("the fitted model's separated group straddles the threshold;",
                       "the planted one is below it")
test_that(title_fixture, {
  m <- load_model(FALSE)
  mu <- as.numeric(m$fitted.values)
  y <- as.numeric(m$y)
  sep <- mu <= 1e-6
  expect_identical(sum(sep), 60L)
  expect_true(all(y[sep] == 0))
  expect_identical(sum(mu <= 1e-8), 13L)
  expect_gt(min(mu[sep]), 1e-9)
  expect_lt(max(mu[sep]), 2e-7)
  expect_gt(min(mu[!sep]), 1e-3)       # four decades of nobody between the group and the rest
  expect_true(isTRUE(m$converged))
  p <- load_model(TRUE)
  mup <- as.numeric(p$fitted.values)
  expect_identical(which(mup <= 1e-8), which(sep))
  expect_equal(mup[sep], rep(1e-12, 60))
  expect_identical(mup[!sep], mu[!sep])
  expect_identical(as.numeric(p$theta), as.numeric(m$theta))
})

test_that("the committed fixture files are what the generator writes", {
  py <- Sys.which("python3")
  skip_if(py == "", "python3 not available")
  tool <- file.path(saige_repo, "harness", "tools", "make_pinned_fixture.py")
  skip_if_not(file.exists(tool))
  status <- suppressWarnings(system2(py, c(tool, "--check"), stdout = FALSE, stderr = FALSE))
  expect_identical(status, 0L)
})

# ---- 3. the planted model ---------------------------------------------------
test_that("planted, BURDEN: guard off the saddlepoint collapses; guard on is the reference (P74)", {
  on <- run_pinned(region_args, dump = TRUE, planted = TRUE, burden = TRUE)
  expect_output_matches_golden(gold("plink_group_pinned_planted_burden.txt"), on$out)
  expect_match(pinned_line(on$log), "^Pinned samples: 60 of 1000 sample", perl = TRUE)
  expect_match(pinned_line(on$log), "left out of the saddlepoint in 2 test\\(s\\)")
  side <- read_side(paste0(on$out, ".spaPinned.txt"))
  expect_identical(side$Kind, c("burden", "burden"))
  expect_identical(side$Group, c("lof", "missense;lof"))
  expect_identical(side$guard, c("applied", "applied"))
  expect_identical(side$n_dropped, c(41L, 43L))
  expect_true(all(side$dropped_mass <= 60e-12 & side$dropped_mass > 0))
  tab_on <- read_saige_table(on$out)
  i <- cells_01(tab_on)
  expect_length(i, 3L)
  # the reference, from the vectors the run saw, two ways
  dumps <- dumps_of(on$dir)
  expect_length(dumps, 2L)
  for (dm in dumps) {
    ref <- reference_p(dm)
    expect_identical(ref$n_material, 10L)
    mc <- monte_carlo_p(dm, ref)
    expect_lt(abs(mc[["p"]] - ref$p), 4 * mc[["se"]])
    expect_true(within(spa_dump_num(dm$h2, "p"), ref$p, 0.15))
    # 6. one quantity two ways: the dump replays to the p the run printed
    expect_identical(sum(dm$d$dropped), if (sum(dm$d$nz) > 100) 43L else 41L)
    expect_true(all(dm$d$gtilde[dm$d$dropped == 1] == 0))
    rp <- replay_spa_dump(dm)
    expect_true(isTRUE(rp$Isconverge))
    expect_equal(rp$pvalue, spa_dump_num(dm$h2, "p"), tolerance = 1e-10)
    expect_true(any(abs(pnum(tab_on$Pvalue_Burden[i]) / rp$pvalue - 1) < 1e-5))
  }
  off <- run_pinned(region_args, env = guard_off, planted = TRUE, burden = TRUE)
  expect_match(pinned_line(off$log), "the guard is OFF")
  expect_false(file.exists(paste0(off$out, ".spaPinned.txt")))
  tab_off <- read_saige_table(off$out)
  expect_identical(dim(tab_off), dim(tab_on))
  # the artefact: 7e-8 and 2e-13 where the reference is 0.02
  expect_true(all(pnum(tab_off$Pvalue_Burden[i]) < 1e-6))
  expect_true(all(pnum(tab_on$Pvalue_Burden[i]) > 0.015 & pnum(tab_on$Pvalue_Burden[i]) < 0.025))
  # nothing else moved: every row but those cells and GENE1's pooled row
  moved <- c(i, which(tab_on$Region == "GENE1" & tab_on$Group == "Cauchy"))
  expect_identical(tab_off[-moved, ], tab_on[-moved, ])
})

test_that("planted, SKAT-O: the saddlepoint R runs itself is guarded too (P75)", {
  on <- run_pinned(region_args, planted = TRUE, burden = FALSE)
  expect_output_matches_golden(gold("plink_group_pinned_planted_skato.txt"), on$out)
  side <- read_side(paste0(on$out, ".spaPinned.txt"))
  expect_setequal(unique(side$Kind), c("region_rescale", "ultra_rare"))
  expect_true(all(side$guard == "applied"))
  expect_true(all(side$Region == "GENE1"))
  # the pseudo-marker's test sees its pinned CARRIERS (the fast kernel); the
  # burden rescale runs the dense kernel on the adjusted burden vector, where
  # every sample has an entry, so all 60 leave it
  expect_true(all(side$n_dropped[side$Kind == "ultra_rare"] == 26L))
  expect_true(all(side$n_dropped[side$Kind == "region_rescale"] == 60L))
  n_tests <- as.integer(sub("^.*left out of the saddlepoint in ([0-9]+) test.*$", "\\1",
                            pinned_line(on$log)))
  expect_identical(n_tests, nrow(side))
  tab_on <- read_saige_table(on$out)
  i <- cells_01(tab_on)
  off <- run_pinned(region_args, env = guard_off, planted = TRUE, burden = FALSE)
  tab_off <- read_saige_table(off$out)
  for (col in c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT")) {
    expect_true(all(pnum(tab_off[[col]][i]) < 1e-10), info = col)
    expect_true(all(pnum(tab_on[[col]][i]) > 0.015 & pnum(tab_on[[col]][i]) < 0.025), info = col)
  }
  moved <- c(i, which(tab_on$Region == "GENE1" & tab_on$Group == "Cauchy"))
  expect_identical(tab_off[-moved, ], tab_on[-moved, ])
})

# ---- 4. the fitted model ----------------------------------------------------
title_fitted <- paste("fitted, BURDEN: the shipped threshold leaves the answer wrong (P72);",
                      "1e-6 does not (P73)")
test_that(title_fitted, {
  shipped <- run_pinned(region_args, dump = TRUE, planted = FALSE, burden = TRUE)
  expect_output_matches_golden(gold("plink_group_pinned_burden.txt"), shipped$out)
  expect_match(pinned_line(shipped$log), "^Pinned samples: 13 of 1000 sample")
  side <- read_side(paste0(shipped$out, ".spaPinned.txt"))
  expect_identical(side$n_dropped, c(10L, 11L))
  expect_identical(side$guard, c("applied", "applied"))
  wide <- run_pinned(region_args, env = eps6, planted = FALSE, burden = TRUE)
  expect_output_matches_golden(gold("plink_group_pinned_burden_eps6.txt"), wide$out)
  expect_match(pinned_line(wide$log), "^Pinned samples: 60 of 1000 sample")
  side6 <- read_side(paste0(wide$out, ".spaPinned.txt"))
  expect_identical(side6$n_dropped, c(41L, 43L))
  off <- run_pinned(region_args, env = guard_off, planted = FALSE, burden = TRUE)
  tabs <- lapply(list(off = off$out, shipped = shipped$out, wide = wide$out), read_saige_table)
  i <- cells_01(tabs$shipped)
  # the two distinct cells have the same score and so the same reference
  refs <- vapply(dumps_of(shipped$dir), function(dm) reference_p(dm)$p, numeric(1))
  expect_length(refs, 2L)
  expect_lt(abs(refs[1] / refs[2] - 1), 0.01)
  ref <- mean(refs)
  p_of <- function(tab) pnum(tab$Pvalue_Burden[i])
  expect_true(all(p_of(tabs$off) / ref > 1.4))
  # CHARACTERIZATION (LEDGER #174, open half): at the shipped threshold 47 of
  # the 60 separated samples stay in the vector and the answer does not come
  # back. Invert this when the threshold is decided.
  expect_true(all(p_of(tabs$shipped) / ref > 1.4))
  expect_true(all(p_of(tabs$shipped) <= p_of(tabs$off)))
  expect_true(all(within(p_of(tabs$wide), ref, 0.15)))
  moved <- c(i, which(tabs$off$Region == "GENE1" & tabs$off$Group == "Cauchy"))
  expect_identical(tabs$off[-moved, ], tabs$shipped[-moved, ])
  expect_identical(tabs$off[-moved, ], tabs$wide[-moved, ])
})

# ---- 5. the budget ----------------------------------------------------------
test_that("over the budget nothing is left out, and the run says so", {
  refused <- run_pinned(region_args, env = c(SAIGE_SPA_PINNED_EPS = "1e-6"),
                        planted = FALSE, burden = TRUE)
  expect_match(pinned_line(refused$log),
               "left out of the saddlepoint in 0 test\\(s\\) and NOT left out in 2 whose")
  side <- read_side(paste0(refused$out, ".spaPinned.txt"))
  expect_identical(side$guard, c("refused", "refused"))
  expect_identical(side$n_dropped, c(0L, 0L))
  expect_true(all(side$dropped_mass > 1e-6))
  off <- run_pinned(region_args, env = guard_off, planted = FALSE, burden = TRUE)
  expect_tables_equivalent(off$out, refused$out, tol = 1e-12)
})

# ---- 7. controls ------------------------------------------------------------
test_that("an unweighted scan does not move when its pinned carriers are left out (P76)", {
  on <- run_pinned(marker_args, planted = TRUE)
  expect_output_matches_golden(gold("plink_marker_pinned_planted.txt"), on$out)
  side <- read_side(paste0(on$out, ".spaPinned.txt"))
  expect_gt(nrow(side), 0)
  expect_true(all(side$Kind == "marker" & side$guard == "applied"))
  off <- run_pinned(marker_args, env = guard_off, planted = TRUE)
  expect_tables_equivalent(off$out, on$out, tol = 1e-6)
})

test_that("a model with no pinned sample reports none and writes no sidecar", {
  out <- tempfile(fileext = ".txt")
  args <- c("--bedFile=./input/genotype_100markers.bed",
            "--bimFile=./input/genotype_100markers.bim",
            "--famFile=./input/genotype_100markers.fam",
            "--AlleleOrder=alt-first", "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
            paste0("--SAIGEOutputFile=", out),
            paste0("--GMMATmodelFile=",
                   file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
            paste0("--varianceRatioFile=",
                   file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt")),
            "--is_output_moreDetails=TRUE")
  st <- run_saige_cli("step2_SPAtests.R", args)
  expect_cli_ok(st)
  log <- readLines(attr(st, "log"), warn = FALSE)
  expect_match(pinned_line(log), "^Pinned samples: 0 \\(no fitted value within 1e-08")
  expect_false(file.exists(paste0(out, ".spaPinned.txt")))
  expect_output_matches_golden(gold("plink_marker_binary.txt"), out)
})

test_that("a setting that does not parse is an error, not a default", {
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", marker_args(out, TRUE),
                      env = c(SAIGE_SPA_PINNED_EPS = "tiny"))
  expect_false(identical(as.integer(st), 0L))
  log <- readLines(attr(st, "log"), warn = FALSE)
  expect_true(any(grepl("SAIGE_SPA_PINNED_EPS='tiny' is not a number in \\[0, 0.5\\)", log)))
})
