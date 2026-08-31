# Internal-consistency invariants over the committed golden outputs. These are
# "computed two ways must reconcile" checks: score-test algebra
# (BETA = Tstat/var, SE = 1/sqrt(var)), allele-count accounting, SPA-vs-normal
# agreement off the SPA path, and the Cauchy combination recomputed in R
# against the C++-written Cauchy rows. Tolerances were measured on the actual
# goldens (limited by 6-sig-fig printing): observed max ~7e-6, asserted 1e-5;
# the ACAT recomputation agrees to ~2e-15, asserted 1e-12.

num <- function(x) suppressWarnings(as.numeric(x))
golden2 <- function(f) read_saige_table(file.path(saige_goldens, "step2", f))

marker_files <- c("plink_marker_binary.txt",
                  "plink_group_fastTest.txt.singleAssoc.txt",
                  "plink_group_cond_binary.txt.singleAssoc.txt")

test_that("score-test algebra: BETA == Tstat/var everywhere", {
  skip_if_no_goldens()
  for (f in c(marker_files, "plink_marker_cond_quantitative.txt")) {
    d <- golden2(f)
    expect_lt(rel_err(num(d$BETA), num(d$Tstat) / num(d$var)), 1e-5)
    if ("BETA_c" %in% colnames(d))
      expect_lt(rel_err(num(d$BETA_c), num(d$Tstat_c) / num(d$var_c)), 1e-5)
  }
})

test_that("score-test algebra: SE == 1/sqrt(var) off the SPA path, and the SPA identity on it", {
  skip_if_no_goldens()
  for (f in marker_files) {
    d <- golden2(f)
    off <- d$Is.SPA == "false"
    expect_gt(sum(off), 0)
    expect_lt(rel_err(num(d$SE[off]), 1 / sqrt(num(d$var[off]))), 1e-5)
    on <- d$Is.SPA == "true"
    if (any(on)) {
      # on SPA rows the reported SE is back-derived from the SPA p-value
      expect_lt(rel_err(num(d$SE[on]),
                        abs(num(d$BETA[on])) / abs(qnorm(num(d$p.value[on]) / 2))), 1e-5)
    }
  }
  dq <- golden2("plink_marker_cond_quantitative.txt")  # no SPA for quantitative
  expect_lt(rel_err(num(dq$SE), 1 / sqrt(num(dq$var))), 1e-5)
})

test_that("p.value equals p.value.NA exactly iff SPA was not applied", {
  skip_if_no_goldens()
  for (f in marker_files) {
    d <- golden2(f)
    off <- d$Is.SPA == "false"; on <- d$Is.SPA == "true"
    expect_true(all(d$p.value[off] == d$p.value.NA[off]))
    if (any(on)) expect_true(all(d$p.value[on] != d$p.value.NA[on]))
  }
})

test_that("allele-count accounting in binary outputs", {
  skip_if_no_goldens()
  for (f in marker_files) {
    d <- golden2(f)
    d <- d[d$CHR != "UR", ]  # ultra-rare collapsed pseudo-markers have no counts
    N <- num(d$N_case) + num(d$N_ctrl)
    expect_lt(rel_err(num(d$AC_Allele2), 2 * num(d$AF_Allele2) * N), 1e-9)  # exact in goldens
    expect_lt(rel_err(num(d$AF_Allele2) * N,
                      num(d$AF_case) * num(d$N_case) + num(d$AF_ctrl) * num(d$N_ctrl)), 1e-5)
  }
})

acat <- function(p) 0.5 - atan(mean(tan((0.5 - p) * pi))) / pi

test_that("Cauchy rows reproduce ACAT over the per-region rows (R recompute vs C++ output)", {
  skip_if_no_goldens()
  # The multi-weight goldens pool over (anno x max_MAF x weight set) rows, so
  # the same one-row-per-region ACAT must hold with K > 1, and it must hold for
  # the additive Pvalue_ACATV/Pvalue_ACATO columns too.
  for (f in c("plink_group_fastTest.txt", "plink_group_cond_binary.txt",
              "plink_group_multiweight.txt", "plink_group_composed.txt")) {
    d <- golden2(f)
    pcols <- intersect(c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
                         "Pvalue_ACATV", "Pvalue_ACATO",
                         "Pvalue_cond", "Pvalue_Burden_cond", "Pvalue_SKAT_cond"),
                       colnames(d))
    for (region in unique(d$Region)) {
      rows <- d[d$Region == region & d$Group != "Cauchy", ]
      cauchy <- d[d$Region == region & d$Group == "Cauchy", ]
      expect_identical(nrow(cauchy), 1L)
      for (col in pcols) {
        expect_lt(rel_err(num(cauchy[[col]]), acat(num(rows[[col]]))), 1e-12)
      }
    }
  }
})

test_that("group-test structural invariants", {
  skip_if_no_goldens()
  d <- golden2("plink_group_fastTest.txt")
  r <- d[d$Group != "Cauchy", ]
  # SKAT-O p lies between (and within 2x of) its component tests
  lo <- pmin(num(r$Pvalue_Burden), num(r$Pvalue_SKAT))
  hi <- pmax(num(r$Pvalue_Burden), num(r$Pvalue_SKAT))
  expect_true(all(num(r$Pvalue) >= lo * (1 - 1e-12) & num(r$Pvalue) <= hi * (1 + 1e-12)))
  expect_true(all(num(r$Pvalue) <= 2 * lo * (1 + 1e-12)))
  # MAC accounting
  expect_equal(num(r$MAC), num(r$MAC_case) + num(r$MAC_control))
  # all-ultra-rare rows degenerate to a single test: string-equal p-values
  allur <- r$Number_rare == "0"
  expect_true(all(r$Pvalue[allur] == r$Pvalue_Burden[allur] &
                  r$Pvalue[allur] == r$Pvalue_SKAT[allur]))
  # Cauchy rows carry NA in the count/effect columns but numeric p-values
  cz <- d[d$Group == "Cauchy", ]
  expect_true(all(cz$max_MAF == "NA" & cz$MAC == "NA"))
  expect_false(anyNA(num(cz$Pvalue)))
})

test_that("conditional group run leaves unconditional columns untouched", {
  skip_if_no_goldens()
  a <- golden2("plink_group_fastTest.txt")
  b <- golden2("plink_group_cond_binary.txt")
  shared <- intersect(colnames(a), colnames(b))
  expect_identical(a[shared], b[shared])
})

# ---- multi-weight goldens (Phase 5 C10) ------------------------------------
# P8 (plink_group_multiweight.txt) is P2's configuration with K = 2 custom sets;
# P8b (plink_group_composed.txt) is P2's group file with two beta pairs plus the
# appended No_Weights set. Both are checked against the SAME P2 golden, so any
# drift in the weight dimension shows up as a structural or arithmetic mismatch
# rather than as an unexplained number.

mw_goldens <- c(plink_group_multiweight.txt = 2L, plink_group_composed.txt = 3L)

test_that("multi-weight goldens: K rows per cell, weight innermost, one Cauchy row", {
  skip_if_no_goldens()
  p2 <- golden2("plink_group_fastTest.txt")
  p2r <- p2[p2$Group != "Cauchy", ]
  for (f in names(mw_goldens)) {
    K <- mw_goldens[[f]]
    d <- golden2(f)
    r <- d[d$Group != "Cauchy", ]
    # the schema columns exist and sit where C7/C8 put them
    expect_identical(colnames(d)[1:5],
                     c("Region", "Group", "max_MAF", "Weight", "Pvalue"))
    expect_true(all(c("Pvalue_ACATV", "Pvalue_ACATO") %in% colnames(d)))
    # exactly K rows per (anno, max_MAF) cell of the single-weight golden
    expect_identical(nrow(r), K * nrow(p2r))
    expect_identical(r$Group, rep(p2r$Group, each = K))
    expect_identical(r$max_MAF, rep(p2r$max_MAF, each = K))
    # weight is the fastest-varying index: the same set order repeats per cell
    sets <- unique(r$Weight)
    expect_identical(length(sets), K)
    expect_identical(r$Weight, rep(sets, nrow(r) / K))
    # one pooled Cauchy row per region, unlabelled by weight
    expect_identical(sum(d$Group == "Cauchy"), length(unique(d$Region)))
    expect_true(all(d$Weight[d$Group == "Cauchy"] == "NA"))
    # MARKER counts are weight-independent -- every set of a cell sees the same
    # markers -- so they repeat across the K rows and match the single-weight
    # golden. (Allele counts are NOT: see the next test.)
    for (col in c("Number_rare", "Number_ultra_rare")) {
      expect_identical(r[[col]], rep(p2r[[col]], each = K))
    }
  }
})

test_that("MAC is set-dependent exactly when the ultra-rare collapse is weighted", {
  skip_if_no_goldens()
  # Design D2: a custom (group-file) set collapses its ultra-rare markers as
  # max(w_r * G), so the pseudo-marker's allele count carries that set's
  # weights and the cell's MAC differs per set. Beta and No_Weights sets
  # collapse UNWEIGHTED, so their MAC is the single-weight golden's.
  p2r <- local({ d <- golden2("plink_group_fastTest.txt"); d[d$Group != "Cauchy", ] })

  custom <- local({ d <- golden2("plink_group_multiweight.txt"); d[d$Group != "Cauchy", ] })
  bySet <- split(num(custom$MAC), custom$Weight)
  expect_identical(length(bySet), 2L)
  expect_false(isTRUE(all.equal(bySet[[1]], bySet[[2]])))
  expect_false(isTRUE(all.equal(bySet[[1]], num(p2r$MAC))))

  composed <- local({ d <- golden2("plink_group_composed.txt"); d[d$Group != "Cauchy", ] })
  bySet2 <- split(num(composed$MAC), composed$Weight)
  expect_identical(length(bySet2), 3L)
  for (s in bySet2) expect_identical(s, num(p2r$MAC))
})

test_that("multi-weight goldens: per-row ACATO is CCT of that row's ACATV/Burden/SKAT", {
  skip_if_no_goldens()
  for (f in names(mw_goldens)) {
    d <- golden2(f)
    r <- d[d$Group != "Cauchy", ]
    for (i in seq_len(nrow(r))) {
      trio <- num(c(r$Pvalue_ACATV[i], r$Pvalue_Burden[i], r$Pvalue_SKAT[i]))
      expect_lt(rel_err(num(r$Pvalue_ACATO[i]), CCT(trio)), 1e-12)
    }
  }
})

test_that("multi-weight goldens: SKAT-O p stays bounded by its components", {
  skip_if_no_goldens()
  for (f in names(mw_goldens)) {
    r <- golden2(f); r <- r[r$Group != "Cauchy", ]
    lo <- pmin(num(r$Pvalue_Burden), num(r$Pvalue_SKAT))
    expect_true(all(num(r$Pvalue) >= lo * (1 - 1e-12)))
    expect_true(all(num(r$Pvalue) <= 2 * lo * (1 + 1e-12)))
    for (col in c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
                  "Pvalue_ACATV", "Pvalue_ACATO")) {
      p <- num(r[[col]])
      expect_false(anyNA(p))
      expect_true(all(p > 0 & p <= 1))
    }
  }
})

test_that("multi-weight BURDEN golden carries the Weight column and matches P8's cells", {
  skip_if_no_goldens()
  b <- golden2("plink_group_multiweight_burden.txt")
  s <- golden2("plink_group_multiweight.txt")
  # written by the C++ BURDEN writer, not the R post-processing, so the schema
  # is asserted independently
  expect_identical(colnames(b)[1:4], c("Region", "Group", "max_MAF", "Weight"))
  expect_false(any(c("Pvalue_ACATV", "Pvalue_ACATO") %in% colnames(b)))
  br <- b[b$Group != "Cauchy", ]; sr <- s[s$Group != "Cauchy", ]
  # same cells and the same weight-innermost ordering as the SKAT-O run
  expect_identical(br$Group, sr$Group)
  expect_identical(br$max_MAF, sr$max_MAF)
  expect_identical(br$Weight, sr$Weight)
  expect_identical(br$MAC, sr$MAC)
})

test_that("composed golden names its sets in the normative order", {
  skip_if_no_goldens()
  d <- golden2("plink_group_composed.txt")
  r <- d[d$Group != "Cauchy", ]
  # beta pairs in the order given on the CLI, then the appended all-ones set
  expect_identical(unique(r$Weight), c("Beta_1_25", "Beta_1_1", "No_Weights"))
})
