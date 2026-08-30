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
  for (f in c("plink_group_fastTest.txt", "plink_group_cond_binary.txt")) {
    d <- golden2(f)
    pcols <- intersect(c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
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
