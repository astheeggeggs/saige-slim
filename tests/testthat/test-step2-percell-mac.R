# LEDGER #118 (the per-set MAC), #126 (the min-group-MAC gate on both region
# paths) and #127 (the pooled row's provenance), end to end through the step-2
# CLI. Tier: standard.
#
# The reconciliation each block runs is named in its comment; nothing here
# compares two different code paths' p-values (the two burden estimators differ
# by 2.76x at MAC 2 -- see LEDGER #126), only cell identity and the accounting
# columns, which are genotype quantities and must agree exactly.

# Aliased at top level, and `saige_extdata` fetched with get(): the idiom
# helper-step2-flexrv.R uses, so that object_usage_linter can see the bindings
# testthat supplies from helper-saige.R at run time.
pc_cli <- run_saige_cli
pc_read <- read_saige_table
pc_ok <- expect_cli_ok
pc_extdata <- function() get("saige_extdata")

PC_G <- c(
  "--bedFile=./input/genotype_100markers.bed",
  "--bimFile=./input/genotype_100markers.bim",
  "--famFile=./input/genotype_100markers.fam",
  "--AlleleOrder=alt-first"
)
PC_MODEL <- c(
  "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
  "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt"
)
PC_SG <- c(
  "--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
  paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
         "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
)

# One run. `groupFile` is an absolute path (the fixtures below are written to
# tempfiles), `extra` the flags that vary.
pc_run <- function(groupFile, extra) {
  out <- file.path(tempfile("step2_pc_"), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  status <- pc_cli("step2_SPAtests.R", c(
    PC_G, PC_MODEL, PC_SG,
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    paste0("--groupFile=", groupFile),
    paste0("--SAIGEOutputFile=", out), extra
  ))
  pc_ok(status)
  list(path = out,
       rows = pc_read(out),
       pooled = paste0(out, ".pooledTests.txt"))
}

pc_cells <- function(r) r$rows[r$rows$Group != "Cauchy", , drop = FALSE]
pc_key <- function(d) {
  paste(d$Region, d$Group, d$max_MAF, if ("Weight" %in% colnames(d)) d$Weight else "", sep = "|")
}

# The bundled two-gene group file, with one custom weight set that is 1 on
# every marker except those in `zero`, where it is exactly 0.
pc_mask_groupfile <- function(zero) {
  g <- withr::with_dir(pc_extdata(), readLines("./input/group_new_snpid.txt"))
  tok <- strsplit(g, " ", fixed = TRUE)
  out <- character(0)
  # var, anno, weight -- in that order. checkGroupFile reads the lines in
  # order and stops with "No anno is specified" if the weight line comes
  # between them.
  for (k in seq(1, length(tok), by = 2)) {
    ids <- tok[[k]][-(1:2)]
    w <- ifelse(ids %in% zero, 0, 1)
    out <- c(out, g[k], g[k + 1],
             paste(c(tok[[k]][1], "weight:MASK", format(w)), collapse = " "))
  }
  f <- tempfile("group_mask_", fileext = ".txt")
  writeLines(out, f)
  f
}

# The same file with `drop`ped markers removed from the var and anno lines --
# the second way of saying "this set does not contain that marker".
pc_drop_groupfile <- function(drop) {
  g <- withr::with_dir(pc_extdata(), readLines("./input/group_new_snpid.txt"))
  tok <- strsplit(g, " ", fixed = TRUE)
  out <- character(0)
  for (k in seq(1, length(tok), by = 2)) {
    ids <- tok[[k]][-(1:2)]
    keep <- !(ids %in% drop)
    var <- c(tok[[k]][1:2], ids[keep])
    anno <- c(tok[[k + 1]][1:2], tok[[k + 1]][-(1:2)][keep])
    out <- c(out, paste(var, collapse = " "), paste(anno, collapse = " "),
             paste(c(tok[[k]][1], "weight:MASK", rep("1", sum(keep))), collapse = " "))
  }
  f <- tempfile("group_drop_", fileext = ".txt")
  writeLines(out, f)
  f
}

PC_BURDEN <- c("--annotation_in_groupTest=lof", "--maxMAF_in_groupTest=0.01",
               "--r.corr=1", "--is_fastTest=FALSE")

test_that("LEDGER #118: weighting a marker 0 equals deleting it, in every column", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # rs4 is a NON-collapsed marker of the GENE1 lof cell (MAC 20 > the default
  # collapse cutoff of 10) -- the arm where the counters were wrong. Weighting
  # it 0 must now agree with removing it from the group file, and the ledger's
  # measurement says the burden STATISTIC already did.
  masked <- pc_run(pc_mask_groupfile("rs4"), PC_BURDEN)
  dropped <- pc_run(pc_drop_groupfile("rs4"), PC_BURDEN)
  kept <- pc_run(pc_mask_groupfile(character(0)), PC_BURDEN)

  m <- pc_cells(masked)
  d <- pc_cells(dropped)
  k <- pc_cells(kept)
  expect_identical(pc_key(m), pc_key(d))
  for (col in c("MAC", "MAC_case", "MAC_control", "Number_rare", "Number_ultra_rare",
                "Pvalue_Burden", "BETA_Burden", "SE_Burden")) {
    expect_identical(m[[col]], d[[col]], label = paste("zeroed vs deleted:", col))
  }
  # POSITIVE CONTROL: the same run with rs4 at weight 1 differs, so the
  # comparison above is not vacuous.
  expect_false(identical(k$MAC, d$MAC))
  expect_false(identical(k$Number_rare, d$Number_rare))
  expect_false(identical(k$Pvalue_Burden, d$Pvalue_Burden))
  # and the drop really removed one rare marker rather than nothing -- from
  # GENE1, which contains rs4, and from GENE1 only, which is the second
  # control: a group file's other genes must not notice.
  expect_identical(d$Region, c("GENE1", "GENE2"))
  expect_identical(as.integer(k$Number_rare) - as.integer(d$Number_rare), c(1L, 0L))
  expect_identical(m$Region[m$MAC != k$MAC], "GENE1")
})

test_that("LEDGER #118: an ultra-rare marker weighted 0 leaves the counter honest", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # rs3 is ultra-rare (MAC 2), so it is inside the collapse. A zero weight
  # already excluded it from the pseudo-marker's VALUE for free -- the paper's
  # weight-before-collapse order -- but `Number_ultra_rare` still counted it.
  masked <- pc_run(pc_mask_groupfile("rs3"), PC_BURDEN)
  dropped <- pc_run(pc_drop_groupfile("rs3"), PC_BURDEN)
  m <- pc_cells(masked)
  d <- pc_cells(dropped)
  expect_identical(m$Number_ultra_rare, d$Number_ultra_rare)
  expect_identical(m$MAC, d$MAC)
  expect_identical(m$Pvalue_Burden, d$Pvalue_Burden)
})

test_that("LEDGER #118: every weight set of a cell reports the plain run's MAC", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # The custom-weight group files have `var`/`anno` lines identical to the
  # unweighted one, so the accounting columns are a property of the CELL and
  # every set of that cell must print them. Before the fix the pseudo-marker's
  # contribution was the WEIGHTED sum, which is what made these columns
  # fractional and set-dependent.
  args <- c("--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
            "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=TRUE")
  plain <- pc_run(withr::with_dir(saige_extdata,
                                  normalizePath("./input/group_new_snpid.txt")), args)
  multi <- pc_run(withr::with_dir(saige_extdata,
                                  normalizePath("./input/group_new_snpid_multiWeights.txt")),
                  args)
  p <- pc_cells(plain)
  m <- pc_cells(multi)
  expect_true("Weight" %in% colnames(m))
  expect_identical(sort(unique(m$Weight)), c("CADD", "REVEL"))
  cols <- c("MAC", "MAC_case", "MAC_control", "Number_rare", "Number_ultra_rare")
  pk <- paste(p$Region, p$Group, p$max_MAF, sep = "|")
  mk <- paste(m$Region, m$Group, m$max_MAF, sep = "|")
  expect_true(all(mk %in% pk))
  idx <- match(mk, pk)
  for (col in cols) {
    expect_identical(m[[col]], p[[col]][idx], label = paste("per-set vs cell:", col))
  }
  # ... and they are integers, which is the point of option (b): a MAC is an
  # allele count, so a gate threshold of 1 or 5 literally means alleles.
  expect_true(all(as.numeric(m$MAC) == round(as.numeric(m$MAC))))
})

# ---------------------------------------------------------------------------
# LEDGER #126: the gate on both paths.
#
# Fixture from the row's own provenance: a region whose single cell has MAC 2.
# --MACCutoff_to_CollapseUltraRare=0 keeps the two singletons as real markers.
# ---------------------------------------------------------------------------

# TWO regions on purpose. GENET is the refused one (MAC 2); GENEU is GENE1's
# lof markers, which no gate in play here refuses. A single-region fixture would
# leave the SKAT-O run with no rows at all and hence a ZERO-BYTE file, which is
# LEDGER #116 and not what this block is about -- and the run would then be
# comparing an unparseable file against a parseable one.
pc_tiny_groupfile <- function() {
  g <- withr::with_dir(pc_extdata(), readLines("./input/group_new_snpid.txt"))
  tok <- strsplit(g, " ", fixed = TRUE)
  gene <- vapply(tok, `[`, "", 1)
  kind <- vapply(tok, `[`, "", 2)
  var <- tok[[which(gene == "GENE1" & kind == "var")]]
  anno <- tok[[which(gene == "GENE1" & kind == "anno")]]
  keep <- anno[-(1:2)] == "lof"
  f <- tempfile("group_tiny_", fileext = ".txt")
  writeLines(c("GENET var rs28 rs32", "GENET anno lof lof",
               paste(c("GENEU", "var", var[-(1:2)][keep]), collapse = " "),
               paste(c("GENEU", "anno", anno[-(1:2)][keep]), collapse = " ")), f)
  f
}

PC_TINY <- c("--annotation_in_groupTest=lof", "--maxMAF_in_groupTest=0.01",
             "--MACCutoff_to_CollapseUltraRare=0", "--is_fastTest=FALSE")

test_that("LEDGER #126: both region paths refuse the same cells at the same gate", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  gf <- pc_tiny_groupfile()
  b5 <- pc_run(gf, c(PC_TINY, "--r.corr=1", "--minGroupMAC_in_BurdenTest=5"))
  s5 <- pc_run(gf, c(PC_TINY, "--r.corr=0", "--minGroupMAC_in_BurdenTest=5"))
  b1 <- pc_run(gf, c(PC_TINY, "--r.corr=1", "--minGroupMAC_in_BurdenTest=1"))
  s1 <- pc_run(gf, c(PC_TINY, "--r.corr=0", "--minGroupMAC_in_BurdenTest=1"))

  # cell identity only -- NEVER the p-values: the two paths do not share a
  # burden estimator and disagree 2.76x at this MAC (LEDGER #126).
  expect_identical(pc_key(pc_cells(b5)), pc_key(pc_cells(s5)))
  expect_identical(pc_key(pc_cells(b1)), pc_key(pc_cells(s1)))
  # GENET (MAC 2) is refused on BOTH paths at 5 and tested on both at 1;
  # GENEU is present throughout, which is what makes the files comparable.
  expect_identical(pc_cells(b5)$Region, "GENEU")
  expect_identical(pc_cells(s5)$Region, "GENEU")
  expect_identical(sort(pc_cells(b1)$Region), c("GENET", "GENEU"))
  expect_identical(sort(pc_cells(s1)$Region), c("GENET", "GENEU"))
  expect_identical(pc_cells(b1)$MAC[pc_cells(b1)$Region == "GENET"], "2")
  expect_identical(pc_cells(s1)$MAC[pc_cells(s1)$Region == "GENET"], "2")
})

test_that("LEDGER #126/#127: a region with nothing tested reports NA, not 0.5", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # CCT of no inputs is the Cauchy CDF at 0, i.e. 0.5 -- a well-formed number
  # that is not a result. The BURDEN emitter writes NA now.
  b5 <- pc_run(pc_tiny_groupfile(), c(PC_TINY, "--r.corr=1",
                                      "--minGroupMAC_in_BurdenTest=5"))
  cauchy <- b5$rows[b5$rows$Group == "Cauchy", , drop = FALSE]
  expect_identical(nrow(cauchy), 2L)   # one per region
  expect_identical(cauchy$Pvalue_Burden[cauchy$Region == "GENET"], "NA")
  # POSITIVE CONTROL: the region that DID test something reports a number.
  expect_false(cauchy$Pvalue_Burden[cauchy$Region == "GENEU"] == "NA")
})

test_that("LEDGER #127: the sidecar names the refused cells and the CCT rule", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  gf <- pc_tiny_groupfile()
  for (rc in c("--r.corr=1", "--r.corr=0")) {
    r5 <- pc_run(gf, c(PC_TINY, rc, "--minGroupMAC_in_BurdenTest=5"))
    expect_true(file.exists(r5$pooled), label = paste("sidecar written", rc))
    tab <- read_saige_table(r5$pooled)
    # only the region with something to disclose is listed; GENEU combined its
    # one cell by the plain Cauchy rule and is absent
    expect_identical(tab$Region, "GENET")
    expect_identical(tab$Cells_refused_minGroupMAC, "1")
    expect_identical(tab$Pvalues_combined, "0")
    expect_identical(tab$Rule, "none")
    expect_identical(tab$Path, if (rc == "--r.corr=1") "BURDEN" else "SKAT-O")
  }
  # NEGATIVE CONTROL: nothing refused, every cell combined by the plain Cauchy
  # rule -> no file at all, and the summary line still prints.
  r1 <- pc_run(gf, c(PC_TINY, "--r.corr=1", "--minGroupMAC_in_BurdenTest=1"))
  expect_false(file.exists(r1$pooled))
})

test_that("LEDGER #126: gate 0 and gate 1 agree byte for byte on both paths", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # Gate 1 is inert BY CONSTRUCTION once the MAC is an allele count: a cell that
  # survives `indexNonZeroVec_arma.n_elem > 0` carries at least one allele. This
  # is what makes the new default safe to ship, so it is asserted rather than
  # assumed -- and it is a claim about ONE path re-run, not two paths compared.
  gf <- withr::with_dir(saige_extdata, normalizePath("./input/group_new_snpid.txt"))
  args <- c("--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
            "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=FALSE")
  for (rc in c("--r.corr=1", "--r.corr=0")) {
    a <- pc_run(gf, c(args, rc, "--minGroupMAC_in_BurdenTest=0"))
    b <- pc_run(gf, c(args, rc, "--minGroupMAC_in_BurdenTest=1"))
    expect_true(files_byte_identical(a$path, b$path), label = paste("gate 0 vs 1", rc))
  }
})
