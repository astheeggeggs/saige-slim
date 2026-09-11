# LEDGER #127, the R seam: which CCT rule a vector of p-values takes, and what
# the sidecar does and does not write. Pure functions and one temp file; no CLI
# run. Tier: fast.
#
# The rule namer has to stay in step with CCT() itself (R/CCT_modified.R) and
# with CCT_cpp (src/CCT.cpp), which apply the same predicates in the same
# order, so the first block checks the naming against the ARITHMETIC rather
# than against a table of expected strings.

test_that("pooledRule names the branch CCT actually takes", {
  expect_identical(pooledRule(numeric(0)), "none")
  expect_identical(pooledRule(c(NA_real_, NA_real_)), "none")
  expect_identical(pooledRule(c(0.2, 0.4)), "Cauchy")
  expect_identical(pooledRule(c(0.2, 1)), "min-p Bonferroni")
  expect_identical(pooledRule(c(0, 0.4)), "zero")
  # a 0 wins over a 1, as CCT's own order does
  expect_identical(pooledRule(c(0, 1)), "zero")
  # NAs are dropped first, exactly as get_CCT_pvalue drops them
  expect_identical(pooledRule(c(NA, 1, 0.3)), "min-p Bonferroni")

  # ... and the naming agrees with what CCT returns (LEDGER #86: a p of exactly
  # 1 turns the combination into min(1, n * min(p)), which is the paper's own
  # prescription and is why the rule has to be disclosed rather than fixed).
  p <- c(0.2, 1)
  expect_identical(CCT(p), min(1, min(p) * length(p)))
  expect_identical(CCT(c(0, 0.4)), 0)
})

test_that("the sidecar is written only when a region has something to disclose", {
  out <- tempfile("pooled_")
  pooledReset()
  # nothing recorded at all: no file, and no summary claim either
  expect_identical(nrow(pooledTable()), 0L)
  expect_silent(pooledSummary(out))
  expect_false(file.exists(paste0(out, ".pooledTests.txt")))

  # a plain region: summary prints, still no file (the negative control has to
  # be visible, so the line is printed either way)
  pooledReset()
  pooledRecord("GENE1", 1L, rows = 4, refused = 0, combined = 4, rule = "Cauchy",
               pvalue = 0.3)
  msg <- capture.output(pooledSummary(out))
  expect_false(file.exists(paste0(out, ".pooledTests.txt")))
  expect_true(any(grepl("plain Cauchy combination", msg)))

  # each of the three disclosure routes on its own writes the file
  for (rec in list(
    list(rows = 4, refused = 2, combined = 4, rule = "Cauchy"),        # gate refused
    list(rows = 4, refused = 0, combined = 3, rule = "Cauchy"),        # an NA dropped
    list(rows = 4, refused = 0, combined = 4, rule = "min-p Bonferroni"),
    list(rows = 0, refused = 0, combined = 0, rule = "none")           # nothing tested
  )) {
    pooledReset()
    do.call(pooledRecord, c(list("GENE1", 1L), rec, list(pvalue = 0.3)))
    path <- paste0(out, ".pooledTests.txt")
    unlink(path)
    invisible(capture.output(pooledSummary(out)))
    expect_true(file.exists(path), label = paste("route", rec$rule, rec$refused, rec$combined))
    tab <- read_saige_table(path)
    expect_identical(nrow(tab), 1L)
    expect_identical(tab$Region, "GENE1")
    unlink(path)
  }
  pooledReset()
})

test_that("the table keeps one row per region, model and pass", {
  pooledReset()
  pooledRecord("GENE1", 1L, 4, 0, 4, "Cauchy", 0.3)
  pooledRecord("GENE1", 2L, 4, 1, 3, "Cauchy", 0.4)
  pooledRecord("GENE2", 1L, 6, 0, 6, "min-p Bonferroni", 1)
  tab <- pooledTable()
  expect_identical(nrow(tab), 3L)
  expect_identical(tab$Model, c(1L, 2L, 1L))
  expect_identical(unique(tab$Path), "SKAT-O")
  expect_identical(colnames(tab),
                   c("Region", "Model", "Path", "Rows", "Pvalues_combined",
                     "Cells_refused_minGroupMAC", "Rule", "Pvalue_pooled"))
  # the C++ half of the table shares that shape, so the two sources bind
  expect_identical(names(getPooledRecordsInCPP()),
                   c("Region", "Model", "Rows", "Pvalues_combined",
                     "Cells_refused_minGroupMAC", "Rule", "Pvalue_pooled"))
  pooledReset()
})
