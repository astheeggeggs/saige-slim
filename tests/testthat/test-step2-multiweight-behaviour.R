# Multi-weight (K > 1) region tests: WEIGHTING BEHAVIOUR AND LEDGER CASES.
#
# Burden linearity under weight scaling, the ultra-rare collapse's effect on
# scale-freeness, the No_Weights set against a plain all-ones run, the forced
# Weight column on an unweighted file, and the characterizations for
# LEDGER #19 (weights.beta reaches C++ but not the R-side ACAT-V) and
# LEDGER #27 (BURDEN + is_fastTest drops every region but the last), plus
# LEDGER #25/#81 (a completed region run validates and is not re-run; FIXED).
#
# Split out of test-step2-multiweight.R; shared scaffolding is in
# helper-multiweight.R. Tier: standard.

test_that("LEDGER #25/#81 FIXED: a completed region run validates its own index and is not re-run", {
  # The writer used to record "nEachChunk = <numberRegionsInChunk>" while that
  # variable was still 0, against a reader expecting 1, so --is_overwrite_output
  # =FALSE failed validation on an untouched index (#25); and it recorded the
  # region counter one past the last region without ever writing the
  # completion line (#81). Both fixed: the index says nEachChunk = 1, its chunk
  # line is the number of regions, the completion line follows, and a re-run
  # with --is_overwrite_output=FALSE says so and changes nothing.
  skip_unless_level("standard")
  out <- file.path(tempfile("step2_mwrestart_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  args <- MW_ARGS("group_new_snpid_withWeights.txt", out)
  expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))
  idx <- readLines(paste0(out, ".index"), warn = FALSE)
  expect_true(any(grepl("^nEachChunk = 1$", idx)))
  expect_identical(tail(idx, 1), "Have completed the analyses of all chunks.")
  expect_identical(tail(idx, 2)[1], "Have completed the analysis of chunk 2")   # two regions
  before <- readLines(out)

  status <- run_saige_cli("step2_SPAtests.R", c(args, "--is_overwrite_output=FALSE"))
  expect_cli_ok(status)
  log <- readLines(attr(status, "log"), warn = FALSE)
  expect_true(any(grepl("the analysis has been completed for the toal 2 chunks", log, fixed = TRUE)))
  expect_false(any(grepl("is not as expected", log)))
  expect_identical(readLines(out), before)
})

test_that("an unweighted file with a forced Weight column labels its single set Beta", {
  res <- implicit_beta_run("1,25")
  expect_true("Weight" %in% colnames(res))
  expect_identical(which(colnames(res) == "Weight"), 4L)   # right after max_MAF
  lab <- unique(res$Weight[res$Group != "Cauchy"])
  expect_identical(lab, "Beta")
  expect_identical(unique(res$Weight[res$Group == "Cauchy"]), "NA")
})

test_that("weights.beta moves the C++ statistics AND the R-side ACAT-V (LEDGER #19 FIXED)", {
  # region_postprocess used to hardcode dbeta(MAF, 1, 25) for the implicit_beta
  # set while C++ honoured g_cfg.weights_beta, so changing --weights.beta
  # shifted every C++-computed column and left Pvalue_ACATV exactly where it
  # was (this test pinned that asymmetry until 2026-09-05). The implicit set now
  # carries the first pair on both sides, so ACAT-V moves with the rest. The
  # cell-for-cell agreement with the explicit-beta path is in
  # test-step2-weightsbeta-single.R.
  a <- implicit_beta_run("1,25")
  b <- implicit_beta_run("1,1")
  expect_identical(colnames(a), colnames(b))

  # R-side weights now follow the pair too
  expect_false(identical(a$Pvalue_ACATV, b$Pvalue_ACATV))
  # computed in C++ from the CLI value, therefore different
  for (cn in c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden", "SE_Burden")) {
    expect_false(identical(a[[cn]], b[[cn]]),
                 label = paste0("weights.beta changes ", cn))
  }
  # allele accounting is weight-independent
  for (cn in c("MAC", "Number_rare", "Number_ultra_rare")) {
    expect_identical(a[[cn]], b[[cn]])
  }
})

# ---- C10 invariants --------------------------------------------------------

test_that("BURDEN + is_fastTest=TRUE keeps every region (LEDGER #27, FIXED)", {
  # Was silent data loss: R signals the fastTest escalation pass by setting
  # flagSparseGRM_cur, but mainRegionInCPP's per-marker variance-ratio
  # selection assigns the SAME flag from MAC, so by the time the BURDEN write
  # decision read it, it reflected the last marker rather than the caller's
  # intent and an escalated region produced no output at all. Only regions
  # taking the cctpval >= 0.1 shortcut survived -- here GENE2 (p = 0.949),
  # while GENE1 (p = 0.0715) escalated and vanished. Main.cpp now snapshots
  # the caller's flag at entry.
  skip_unless_level("standard")
  for (gf in c("group_new_snpid.txt", "group_new_snpid_multiWeights.txt")) {
    out <- file.path(tempfile("step2_burdenfast_"), "region.txt")
    dir.create(dirname(out), recursive = TRUE)
    args <- MW_ARGS(gf, out)
    args <- c(args[!grepl("^--is_fastTest=", args)], "--r.corr=1", "--is_fastTest=TRUE")
    expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))
    d <- read_saige_table(out)
    expect_identical(sort(unique(d$Region)), c("GENE1", "GENE2"))
  }
})

test_that("BURDEN fastTest and non-fastTest agree row for row (LEDGER #27)", {
  # The point of fastTest is to escalate when the cheap pass is not good
  # enough, so its output must match the accurate path. Before the fix this
  # could not even be asked -- the fast run was missing a region.
  skip_unless_level("standard")
  runs <- lapply(c("TRUE", "FALSE"), function(ft) {
    out <- file.path(tempfile(paste0("step2_burdencmp_", ft, "_")), "region.txt")
    dir.create(dirname(out), recursive = TRUE)
    args <- MW_ARGS("group_new_snpid.txt", out)
    args <- c(args[!grepl("^--is_fastTest=", args)], "--r.corr=1",
              paste0("--is_fastTest=", ft))
    expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))
    read_saige_table(out)
  })
  fast <- runs[[1]][order(runs[[1]]$Region, runs[[1]]$Group, runs[[1]]$max_MAF), ]
  slow <- runs[[2]][order(runs[[2]]$Region, runs[[2]]$Group, runs[[2]]$max_MAF), ]
  expect_identical(nrow(fast), nrow(slow))
  expect_identical(fast$Region, slow$Region)
  expect_identical(fast$Group, slow$Group)
  # statistics at ~1e-12 relative, never bit-identity across two code paths
  for (cc in c("Pvalue_Burden", "BETA_Burden", "SE_Burden")) {
    a <- suppressWarnings(as.numeric(fast[[cc]]))
    b <- suppressWarnings(as.numeric(slow[[cc]]))
    k <- !is.na(a) & !is.na(b)
    expect_equal(a[k], b[k], tolerance = 1e-12, label = cc)
  }
})

test_that("BURDEN + is_fastTest=FALSE keeps every region", {
  skip_unless_level("standard")
  out <- file.path(tempfile("step2_burdenslow_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  args <- MW_ARGS("group_new_snpid_multiWeights.txt", out)
  args <- c(args[!grepl("^--is_fastTest=", args)], "--r.corr=1", "--is_fastTest=FALSE")
  expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))
  d <- read_saige_table(out)
  expect_identical(sort(unique(d$Region)), c("GENE1", "GENE2"))
})

test_that("burden linearity: scaling every weight by 2 is exact without the UR collapse", {
  # A uniform rescaling of user-supplied weights is a reparameterisation of the
  # burden test: the score scales by c, its variance by c^2, so the p-value is
  # unchanged and BETA/SE both scale by 1/c. Checked with the ultra-rare
  # collapse DISABLED, where the identity is exact to machine precision.
  skip_unless_level("standard")
  doubled <- local({
    src <- file.path(saige_extdata, "input", "group_new_snpid_multiWeights.txt")
    lines <- readLines(src)
    out <- tempfile(fileext = ".txt")
    writeLines(vapply(lines, function(l) {
      p <- strsplit(l, "[ \t]+")[[1]]
      if (length(p) >= 3 && grepl("^weight", p[2])) {
        p[-(1:2)] <- format(as.numeric(p[-(1:2)]) * 2, trim = TRUE, scientific = FALSE)
        paste(p, collapse = " ")
      } else l
    }, "", USE.NAMES = FALSE), out)
    out
  })
  run_one <- function(groupFile) {
    out <- file.path(tempfile("step2_lin_"), "region.txt")
    dir.create(dirname(out), recursive = TRUE)
    args <- MW_ARGS("group_new_snpid_multiWeights.txt", out)
    args[grep("^--groupFile=", args)] <- paste0("--groupFile=", groupFile)
    args <- c(args, "--MACCutoff_to_CollapseUltraRare=0")
    expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))
    read_saige_table(out)
  }
  a <- run_one(file.path(saige_extdata, "input", "group_new_snpid_multiWeights.txt"))
  b <- run_one(doubled)
  ra <- a[a$Group != "Cauchy", ]
  rb <- b[b$Group != "Cauchy", ]
  expect_identical(nrow(ra), nrow(rb))
  expect_true(all(pnum(ra$Number_ultra_rare) == 0))    # collapse really is off
  # scale-invariant: the p-values are unchanged, to the printed precision
  # toleranced rather than string-identical: the two runs feed different
  # intermediate magnitudes (every weight doubled) into the same arithmetic, so
  # the last ULP need not agree even though the p-value is invariant
  expect_lt(rel_err(pnum(ra$Pvalue_Burden), pnum(rb$Pvalue_Burden)), 1e-12)
  expect_lt(rel_err(pnum(ra$Pvalue_SKAT), pnum(rb$Pvalue_SKAT)), 1e-12)
  # equivariant: BETA and SE both halve
  expect_lt(max(abs(pnum(rb$BETA_Burden) / pnum(ra$BETA_Burden) - 0.5)), 1e-9)
  expect_lt(max(abs(pnum(rb$SE_Burden) / pnum(ra$SE_Burden) - 0.5)), 1e-9)
})

test_that("the ultra-rare collapse makes custom weights only approximately scale-free", {
  # With the collapse ON the same doubling stops being exact: the collapsed
  # pseudo-genotype is max(w_r*G), so its apparent allele count doubles too
  # (AC 38.40 -> 76.80 in the singleAssoc output) while its weight in
  # AnnoWeights stays 1, and the MAC-dependent machinery downstream therefore
  # sees a different marker. So a uniform rescaling of a group file -- which a
  # user would reasonably expect to be a no-op -- shifts the answer slightly.
  # Bounded, not asserted equal: a real weight-handling regression would move
  # these by orders of magnitude, while this residual is ~1e-2.
  skip_unless_level("standard")
  doubled <- local({
    lines <- readLines(file.path(saige_extdata, "input",
                                 "group_new_snpid_multiWeights.txt"))
    out <- tempfile(fileext = ".txt")
    writeLines(vapply(lines, function(l) {
      p <- strsplit(l, "[ \t]+")[[1]]
      if (length(p) >= 3 && grepl("^weight", p[2])) {
        p[-(1:2)] <- format(as.numeric(p[-(1:2)]) * 2, trim = TRUE, scientific = FALSE)
        paste(p, collapse = " ")
      } else l
    }, "", USE.NAMES = FALSE), out)
    out
  })
  out <- file.path(tempfile("step2_lin_on_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  args <- MW_ARGS("group_new_snpid_multiWeights.txt", out)
  args[grep("^--groupFile=", args)] <- paste0("--groupFile=", doubled)
  expect_cli_ok(run_saige_cli("step2_SPAtests.R", args))

  a <- multi_run()$region                      # collapse on, original weights
  b <- read_saige_table(out)                   # collapse on, doubled weights
  ra <- a[a$Group != "Cauchy", ]
  rb <- b[b$Group != "Cauchy", ]
  expect_true(all(pnum(ra$Number_ultra_rare) > 0))   # collapse really is active
  rel <- function(x, y) max(abs(x - y) / (abs(x) + abs(y) + 1e-300))
  # not exact any more ...
  expect_false(identical(ra$Pvalue_Burden, rb$Pvalue_Burden))
  # ... but still close, and BETA/SE still halve to within a couple of percent
  expect_lt(rel(pnum(ra$Pvalue_Burden), pnum(rb$Pvalue_Burden)), 0.05)
  expect_lt(max(abs(pnum(rb$BETA_Burden) / pnum(ra$BETA_Burden) - 0.5)), 0.05)
  expect_lt(max(abs(pnum(rb$SE_Burden) / pnum(ra$SE_Burden) - 0.5)), 0.05)
})

test_that("No_Weights rows reproduce a plain all-ones run", {
  # The appended all-ones set must be the same analysis as running the
  # unweighted group file in no-weight mode. This is what pins design D2's
  # "noweight = all ones with an unweighted UR collapse" against the C++ side.
  #
  # Counts and labels are compared exactly; the statistics are compared at a
  # tolerance because the two answers arrive by different routes -- the
  # No_Weights column is the third set of a K = 3 run, the plain run is K = 1 --
  # so they are algebraically identical but not bit-identical. Observed
  # disagreement is 1 ULP on 1 of 12 rows (0.540776608214742 vs ...741), and
  # whether it shows up in the printed digits is platform-dependent: it does
  # not on the golden-generating macOS host, it does on the Linux CI runner.
  # Asserting string identity here is over-claiming; 1e-12 still catches any
  # real divergence between the two paths by many orders of magnitude.
  skip_unless_level("standard")
  plain_out <- file.path(tempfile("step2_nwplain_"), "region.txt")
  dir.create(dirname(plain_out), recursive = TRUE)
  expect_cli_ok(run_saige_cli(
    "step2_SPAtests.R",
    c(MW_ARGS("group_new_snpid.txt", plain_out), "--is_no_weight_in_groupTest=TRUE")
  ))
  plain <- read_saige_table(plain_out)

  d <- noweights_run()$region                  # weighted file + No_Weights set
  nw <- d[d$Weight == "No_Weights" & d$Group != "Cauchy", ]
  pr <- plain[plain$Group != "Cauchy", ]
  expect_identical(nrow(nw), nrow(pr))
  expect_identical(nw$Group, pr$Group)
  expect_identical(nw$max_MAF, pr$max_MAF)
  # exact: labels, marker counts and allele counts are integers or short decimals
  for (col in c("MAC", "MAC_case", "MAC_control", "Number_rare",
                "Number_ultra_rare")) {
    expect_identical(nw[[col]], pr[[col]], label = paste("No_Weights ==", col))
  }
  # toleranced: the fitted statistics
  for (col in c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden",
                "SE_Burden")) {
    expect_lt(rel_err(pnum(nw[[col]]), pnum(pr[[col]])), 1e-12)
  }
})
