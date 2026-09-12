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
              "plink_group_multiweight.txt", "plink_group_composed.txt",
              "plink_group_cond_multiweight.txt")) {
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
  # MAC accounting. Holds here because no cell of this golden sits above AF
  # 0.5; it is NOT a general identity -- the group MAC is folded and the
  # case/control split is not, so the two diverge on any cell above 0.5
  # (LEDGER #133, exhibited by P62/P63 and pinned in test-step2-region-flip.R).
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
  # the same invariant at K = 2 (C11, design D7): P10 is P8 plus --condition,
  # and the per-set conditional adjustment must not move any unconditional
  # column (Pvalue_ACATV/Pvalue_ACATO included)
  a <- golden2("plink_group_multiweight.txt")
  b <- golden2("plink_group_cond_multiweight.txt")
  shared <- intersect(colnames(a), colnames(b))
  expect_identical(a[shared], b[shared])
  # positive control on the conditional side: the adjustment is per set, so
  # the two sets' conditional p-values differ from each other and from the
  # unconditional ones
  r <- b[b$Group != "Cauchy", ]
  expect_false(identical(r$Pvalue_cond[r$Weight == "CADD"],
                         r$Pvalue_cond[r$Weight == "REVEL"]))
  expect_false(identical(r$Pvalue_cond, r$Pvalue))
})

# ---- multi-weight goldens (Phase 5 C10) ------------------------------------
# P8 (plink_group_multiweight.txt) is P2's configuration with K = 2 custom sets;
# P8b (plink_group_composed.txt) is P2's group file with two beta pairs plus the
# appended No_Weights set. Both are checked against the SAME P2 golden, so any
# drift in the weight dimension shows up as a structural or arithmetic mismatch
# rather than as an unexplained number.

mw_goldens <- c(plink_group_multiweight.txt = 2L, plink_group_composed.txt = 3L,
                plink_group_cond_multiweight.txt = 2L)

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

test_that("MAC is a property of the CELL, whatever the weight set (LEDGER #118)", {
  skip_if_no_goldens()
  # INVERTED 2026-09-11 with the LEDGER #118 fix. It used to read "MAC is
  # set-dependent exactly when the ultra-rare collapse is weighted": a custom
  # (group-file) set collapses its ultra-rare markers as max(w_r * G), and the
  # pseudo-marker's WEIGHTED allele count went into the MAC column, so the same
  # cell printed a different MAC for each set (38.4 for CADD against 3.34 for
  # REVEL) and neither matched the unweighted run. Under option (b) the column
  # is an allele count throughout -- the pseudo-marker contributes its
  # UNWEIGHTED collapse -- so every set of a cell prints the cell's own numbers.
  # The group files make this exact: group_new_snpid_multiWeights.txt has
  # `var`/`anno` lines byte-identical to group_new_snpid.txt, so this is one
  # quantity computed two ways and equality is the only correct answer.
  p2r <- local({ d <- golden2("plink_group_fastTest.txt"); d[d$Group != "Cauchy", ] })

  custom <- local({ d <- golden2("plink_group_multiweight.txt"); d[d$Group != "Cauchy", ] })
  bySet <- split(num(custom$MAC), custom$Weight)
  expect_identical(length(bySet), 2L)
  expect_identical(bySet[[1]], bySet[[2]])
  for (s in bySet) expect_identical(unname(s), num(p2r$MAC))
  # ... and it is an allele count, not a weighted sum: every value an integer
  expect_identical(num(custom$MAC), round(num(custom$MAC)))

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

test_that("ACAT-V recomputes from the single-variant p-values and the fixture weights", {
  skip_if_no_goldens()
  # The strongest available check on Pvalue_ACATV: rebuild each cell's ACAT-V
  # from first principles -- the group file's per-set weights, the cell's
  # annotation/max_MAF membership, the per-marker p-values in the singleAssoc
  # golden, and weight 1 for the collapsed ultra-rare pseudo-marker (design
  # D2's custom-set rule, since C++ has already baked w_r into the collapse).
  # The per-row ACATO check above consumes ACATV as an input, so without this
  # nothing pins ACATV's own value.
  reg <- golden2("plink_group_multiweight.txt")
  sa  <- golden2("plink_group_multiweight.txt.singleAssoc.txt")
  gl <- strsplit(readLines(file.path(saige_extdata, "input",
                                     "group_new_snpid_multiWeights.txt")), "[ \t]+")
  grouprow <- function(gene, pattern) {
    i <- which(vapply(gl, function(p) p[1] == gene && grepl(pattern, p[2]), TRUE))[1]
    gl[[i]][-(1:2)]
  }

  r <- reg[reg$Group != "Cauchy", ]
  expect_gt(nrow(r), 0)
  for (i in seq_len(nrow(r))) {
    cell_id <- paste(r$Region[i], r$Group[i], r$max_MAF[i], r$Weight[i], sep = ":")
    p_ur <- num(sa$p.value[sa$MarkerID == cell_id])
    expect_identical(length(p_ur), 1L)   # one UR pseudo-marker per cell

    if (r$Number_rare[i] == "0") {
      # nothing but the pseudo-marker in the cell: ACAT-V is its p-value
      expect_identical(num(r$Pvalue_ACATV[i]), p_ur)
      next
    }
    vars  <- grouprow(r$Region[i], "^var$")
    annos <- grouprow(r$Region[i], "^anno$")
    w     <- as.numeric(grouprow(r$Region[i], paste0("weight.", r$Weight[i])))
    real  <- sa[sa$CHR != "UR", ]
    idx   <- match(vars, real$MarkerID)
    keep  <- annos %in% strsplit(r$Group[i], ";")[[1]] &
             !is.na(idx) & num(real$AF_Allele2[idx]) <= num(r$max_MAF[i])
    expect_identical(sum(keep), as.integer(num(r$Number_rare[i])))
    expect_lt(rel_err(num(r$Pvalue_ACATV[i]),
                      CCT(c(num(real$p.value[idx[keep]]), p_ur), c(w[keep], 1))),
              1e-12)
  }
})

test_that("useSparseGRMforVarRatio changes only the variance-ratio file", {
  # binary_sparseGRM_vr is binary_sparseGRM_a's step-1 command plus
  # --useSparseGRMforVarRatio=TRUE. The flag adds `sparse` rows to the variance
  # ratio file and touches nothing else: the fitted model is byte-identical and
  # the null / null_noXadj rows are unchanged. That is what makes P27-P29
  # diffable against P1 -- if this fails, those diffs stop isolating the step-2
  # branch and start including a different model.
  skip_if_no_goldens()
  s1 <- function(f) file.path(saige_goldens, "step1", f)
  skip_if_not(file.exists(s1("binary_sparseGRM_vr.rda")), "vr step-1 golden absent")
  expect_identical(tools::md5sum(s1("binary_sparseGRM_vr.rda"))[[1]],
                   tools::md5sum(s1("binary_sparseGRM_a.rda"))[[1]])

  vr <- read.table(s1("binary_sparseGRM_vr.varianceRatio.txt"))
  a <- read.table(s1("binary_sparseGRM_a.varianceRatio.txt"))
  expect_identical(sort(unique(vr$V2)), c("null", "null_noXadj", "sparse"))
  # the non-sparse rows are the same numbers, types and categories, in order
  keep <- vr$V2 != "sparse"
  expect_identical(vr$V1[keep], a$V1)
  expect_identical(vr$V2[keep], a$V2)
  expect_identical(vr$V3[keep], a$V3)
  # sparse rows: one per MAC category, and near 1 (full-GRM numerator over
  # sparse-GRM denominator on a fixture whose relatedness is nearly all captured)
  expect_identical(sum(vr$V2 == "sparse"), sum(a$V2 == "null"))
  expect_true(all(abs(vr$V1[vr$V2 == "sparse"] - 1) < 1e-3))
})

test_that("the same holds for the QUANTITATIVE variance-ratio pair", {
  # quant_sparseGRM_vr is quant_sparseGRM's step-1 command plus
  # --useSparseGRMforVarRatio=TRUE, and it is what let the quantitative
  # sparse-GRM step-2 path get a golden at all (P38-P41, 2026-09-08): a region
  # test needs `sparse` rows and only that flag writes them. The same invariant
  # as the binary pair above, for the same reason -- it is what makes P38/P39
  # diffable against P37/P36.
  #
  # It is also the only golden containing a quantitative `sparse` row, so it is
  # where LEDGER #97 is visible: #65 scales the null / null_noXadj ratios by
  # 1/s^2 on the way out of a quantitative fit and leaves `sparse` alone,
  # because that ratio is scale-free.
  skip_if_no_goldens()
  s1 <- function(f) file.path(saige_goldens, "step1", f)
  skip_if_not(file.exists(s1("quant_sparseGRM_vr.rda")), "quant vr step-1 golden absent")
  expect_identical(tools::md5sum(s1("quant_sparseGRM_vr.rda"))[[1]],
                   tools::md5sum(s1("quant_sparseGRM.rda"))[[1]])

  vr <- read.table(s1("quant_sparseGRM_vr.varianceRatio.txt"))
  a <- read.table(s1("quant_sparseGRM.varianceRatio.txt"))
  expect_identical(sort(unique(vr$V2)), c("null", "null_noXadj", "sparse"))
  keep <- vr$V2 != "sparse"
  expect_identical(vr$V1[keep], a$V1)
  expect_identical(vr$V2[keep], a$V2)
  expect_identical(vr$V3[keep], a$V3)
  expect_identical(sum(vr$V2 == "sparse"), sum(a$V2 == "null"))
  # Near 1, but NOT as near as the binary pair's 1e-3: measured 0.9969 and
  # 0.9976, i.e. ~3.2e-3 out. That gap is the same physics that makes this
  # golden worth having. The sparse ratio is a full-GRM numerator over a
  # sparse-GRM denominator, so it departs from 1 in proportion to the GRM's
  # share of Sigma = tau0/W + tau1*GRM. Here W = 1 and tau = (0.216, 0.465), so
  # the GRM term is 2.2x the diagonal one; on the binary fit tau0/W >= 4 against
  # tau1 = 0.33, so it is at most 0.08 -- about 26x less. The same ordering is
  # why LEDGER #106's one-triangle Sigma was first-order on a quantitative trait
  # and ~0.5% on a binary one.
  expect_true(all(abs(vr$V1[vr$V2 == "sparse"] - 1) < 1e-2))
  expect_true(all(abs(vr$V1[vr$V2 == "sparse"] - 1) > 1e-3))
})

# ---- FlexRV goldens (F4, 2026-09-09; P58-P61) --------------------------------
# One cell (lof;missense at max MAF 0.01) and the default 192 transform sets on
# group_new_snpid_flexrv.txt. The fixture arithmetic is in
# helper-step2-flexrv.R: at the model's 1,000 samples no in-cell marker has
# MAC > 20, so MAC25..MAC500 are duplicates of unw; at the default collapse
# cutoff MAC1/3/5/10 select only collapsed markers.

flexrv_goldens <- c(plink_group_flexrv_burden.txt = 10L, plink_group_flexrv_burden_quant.txt = 0L)

test_that("FlexRV goldens: 192 sets per gene in the normative order, one Cauchy row, p_FlexRV = CCT", {
  skip_if_no_goldens()
  names192 <- flexrvTransformTable(flexrvDefaultScoreSpec, flexrvDefaultMAFSpec, "REVEL",
                                   maxMAF = 0.01)$names
  for (f in names(flexrv_goldens)) {
    d <- golden2(f)
    expect_identical(colnames(d)[1:5], c("Region", "Group", "max_MAF", "Weight", "Pvalue_Burden"))
    r <- d[d$Group != "Cauchy", ]
    cau <- d[d$Group == "Cauchy", ]
    expect_identical(cau$Region, c("GENE1", "GENE2"))
    for (g in c("GENE1", "GENE2")) {
      rg <- r[r$Region == g, ]
      expect_identical(rg$Weight, names192, label = paste(f, g))
      expect_true(all(rg$Group == "lof;missense" & rg$max_MAF == "0.01"))
      p <- num(rg$Pvalue_Burden)
      expect_true(all(is.finite(p) & p > 0 & p <= 1))
      # the pooled row is the CCT over the 192 rows, from their printed digits
      expect_lt(rel_err(num(cau$Pvalue_Burden[cau$Region == g]), CCT(p)), 1e-5)
      # the ultra-rare counts say which run this is: collapse on (10) or off (0)
      expect_true(all(num(rg$Number_ultra_rare) >= (if (flexrv_goldens[[f]] > 0) 1 else 0)))
      if (flexrv_goldens[[f]] == 0) expect_true(all(num(rg$Number_ultra_rare) == 0))
    }
  }
})

test_that("FlexRV goldens: a MAC threshold above every MAC is unw, and the axis is otherwise live", {
  skip_if_no_goldens()
  strip <- function(d) {
    d$Weight <- NULL
    rownames(d) <- NULL
    d
  }
  # measured 2026-09-09 when the goldens were made: the number of DISTINCT
  # (Pvalue_Burden, BETA_Burden, SE_Burden) triples per gene among the 192
  # rows. Not 192, because MAC25..MAC500 repeat unw (6 of 12 MAF transforms
  # coincide) and, with the collapse on, the four low MAC sets differ only
  # through the pseudo-marker. A change here is a change in which sets the
  # fixture can tell apart, which is exactly what a coverage claim rests on.
  # GENE2's count moved 102 -> 104 with the LEDGER #118 fix: two sets whose
  # per-set MAC now lands in a different variance-ratio band from the cell's
  # became distinguishable. It is a COVERAGE measurement of this fixture, so it
  # is re-measured when the numbers legitimately move, never predicted.
  #
  # LEDGER #129 moved the collapse-ON counts 100 -> 111 and 104 -> 111, and
  # WHERE they landed is the point: 111 is what the collapse-OFF run has always
  # had. The dosage filter had been cleaning several sets' weighted
  # pseudo-markers down to one shared support -- REVEL_cbrt0.4..0.9_MAC1 all
  # printed the same 1.453839E-02 -- so the fixture's power to tell sets apart
  # depended on --MACCutoff_to_CollapseUltraRare. It no longer does, and the two
  # rows of this list agreeing is that statement. Re-measured on the fixed
  # build, not carried over from the screen that predicted it.
  distinct <- list(plink_group_flexrv_burden.txt = c(GENE1 = 111L, GENE2 = 111L),
                   plink_group_flexrv_burden_quant.txt = c(GENE1 = 111L, GENE2 = 111L))
  for (f in names(flexrv_goldens)) {
    d <- golden2(f)
    r <- d[d$Group != "Cauchy", ]
    for (g in c("GENE1", "GENE2")) {
      rg <- r[r$Region == g, ]
      for (sc in c("unw", "lof", "cbrt0.2", "sbeta20")) {
        base <- strip(rg[rg$Weight == paste0("REVEL_", sc, "_unw"), ])
        for (t in c(25, 50, 100, 200, 500)) {
          expect_identical(strip(rg[rg$Weight == paste0("REVEL_", sc, "_MAC", t), ]), base,
                           label = paste(f, g, sc, t))
        }
        expect_false(identical(strip(rg[rg$Weight == paste0("REVEL_", sc, "_MAC10"), ]), base))
        expect_false(identical(strip(rg[rg$Weight == paste0("REVEL_", sc, "_mafbeta2@0.01"), ]),
                               base))
      }
      key <- paste(rg$Pvalue_Burden, rg$BETA_Burden, rg$SE_Burden)
      expect_identical(length(unique(key)), distinct[[f]][[g]], label = paste(f, g))
    }
  }
})

test_that("FlexRV goldens: the pgen mirror is byte-identical to its PLINK twin", {
  skip_if_no_goldens()
  # the MAC ladder and the mafbeta transforms read the IN-SAMPLE MAC/MAF, which
  # must not depend on the reader
  expect_true(files_byte_identical(file.path(saige_goldens, "step2", "plink_group_flexrv_burden.txt"),
                                   file.path(saige_goldens, "step2", "pgen_group_flexrv_burden.txt")))
})

test_that("FlexRV SKAT-O golden: 4 sets per gene, multi-weight schema, no kernel failure", {
  skip_if_no_goldens()
  d <- golden2("plink_group_flexrv_skato.txt")
  expect_identical(colnames(d)[1:5], c("Region", "Group", "max_MAF", "Weight", "Pvalue"))
  expect_true(all(c("Pvalue_ACATV", "Pvalue_ACATO") %in% colnames(d)))
  r <- d[d$Group != "Cauchy", ]
  expect_identical(r$Weight, rep(c("REVEL_unw_unw", "REVEL_unw_mafbeta2@0.02", "REVEL_sbeta2_unw",
                                   "REVEL_sbeta2_mafbeta2@0.02"), 2))
  expect_false(any(r$Pvalue == "NA"))
  for (i in seq_len(nrow(r))) {
    trio <- num(c(r$Pvalue_ACATV[i], r$Pvalue_Burden[i], r$Pvalue_SKAT[i]))
    expect_lt(rel_err(num(r$Pvalue_ACATO[i]), CCT(trio)), 1e-12)
  }
  expect_identical(sum(d$Group == "Cauchy"), 2L)
})
