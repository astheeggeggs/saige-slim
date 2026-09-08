# The QUANTITATIVE region path (goldens P36/P37), which had no coverage at all:
# every other group golden is binary, so the `traitType == "quantitative"` arms
# of the C++ region loop and of region_postprocess were gated by nothing.
#
# P36 is the run LEDGER #72 needs -- `--r.corr=1` puts a quantitative trait
# through the BURDEN Cauchy aggregation whose `std::stod` is bare (it cannot
# make that parse throw here; what it adds is a trait type with no SPA rescue
# arm on the line). P37 is the default SKAT-O arm, and it is what surfaced
# LEDGER #86.
#
# Both models use quant_sparseGRM WITHOUT --sparseGRMFile: it is a
# --useSparseGRMtoFitNULL fit, so its variance-ratio file carries no `sparse`
# rows and step 2 must not be given the sparse GRM for it.

goldens_step2_dir <- file.path(saige_goldens, "step2")

# Top level, not inside a function: `saige_goldens` comes from helper-saige.R,
# which testthat sources first.
quant_step1 <- file.path(saige_goldens, "step1", "quant_sparseGRM")
quant_model <- c(paste0("--GMMATmodelFile=", quant_step1, ".rda"),
                 paste0("--varianceRatioFile=", quant_step1, ".varianceRatio.txt"))

golden_step2 <- function(f) file.path(goldens_step2_dir, f)

quant_group_args <- function(out) {
  c("--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first",
    paste0("--SAIGEOutputFile=", out),
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    quant_model,
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01")
}

test_that("P36: quantitative BURDEN group test reproduces its golden", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- tempfile(fileext = ".txt")
  expect_cli_ok(run_saige_cli("step2_SPAtests.R",
                              c(quant_group_args(out), "--r.corr=1", "--is_fastTest=FALSE")))
  gold <- golden_step2("plink_group_burden_quant.txt")
  expect_output_matches_golden(gold, out)

  # It is the BURDEN schema (no SKAT-O columns) and it does carry the Cauchy
  # rows that the bare std::stod at src/Main.cpp:2793 feeds -- LEDGER #72.
  d <- read_saige_table(out)
  expect_identical(colnames(d)[1:4], c("Region", "Group", "max_MAF", "Pvalue_Burden"))
  expect_false("Pvalue_SKAT" %in% colnames(d))
  expect_identical(sum(d$Group == "Cauchy"), 2L)
  expect_false(anyNA(suppressWarnings(as.numeric(d$Pvalue_Burden[d$Group == "Cauchy"]))))
})

test_that("P37: quantitative SKAT-O group test reproduces its golden", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- tempfile(fileext = ".txt")
  expect_cli_ok(run_saige_cli("step2_SPAtests.R",
                              c(quant_group_args(out), "--is_fastTest=FALSE")))
  gold <- golden_step2("plink_group_nofast_quant.txt")
  expect_output_matches_golden(gold, out)
  expect_output_matches_golden(paste0(gold, ".singleAssoc.txt"),
                               paste0(out, ".singleAssoc.txt"))
})

# --- P38-P41: the quantitative SPARSE-GRM path (added 2026-09-08) -----------
# P36/P37 above deliberately do NOT pass --sparseGRMFile, because
# quant_sparseGRM's variance-ratio file has no `sparse` rows. The consequence,
# measured over all 38 step-2 goldens before these were added: 23 of the 33
# binary runs passed --sparseGRMFile and 0 of the 5 quantitative ones did, so
# NOTHING put a quantitative model through setSparseSigma_new or
# getPCG1ofSigmaAndGtilde.
#
# That is the regime where LEDGER #106 mattered most, by its own WHAT cell: the
# one-triangle Sigma is "small here only because 1/W ~ 11 dominates
# tau2 x kinship ~ 0.16; on a QUANTITATIVE trait (diagonal tau1 ~ 1) ... the
# missing triangle is first-order". Measured on P39 against the pre-#106 build:
# BETA_Burden moved by up to 64%, Pvalue_Burden by up to 12%, against ~0.5% on
# the binary goldens. The #106 fix moved 43 binary goldens and could not move a
# single quantitative byte, because there were none to move.
#
# quant_sparseGRM_vr is quant_sparseGRM's step-1 command plus
# --useSparseGRMforVarRatio=TRUE, which is the only way to get `sparse` rows;
# test-invariants-goldens.R pins that the .rda is byte-identical to
# quant_sparseGRM.rda, so a diff of P38/P39 against P37/P36 isolates the step-2
# branch and not the model.

quant_vr_step1 <- file.path(saige_goldens, "step1", "quant_sparseGRM_vr")
quant_vr_model <- c(paste0("--GMMATmodelFile=", quant_vr_step1, ".rda"),
                    paste0("--varianceRatioFile=", quant_vr_step1, ".varianceRatio.txt"))
quant_sparse_grm <- c(
  "--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx",
  paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed",
         ".sparseGRM.mtx.sampleIDs.txt"))

quant_sparse_group_args <- function(out, geno) {
  c(geno,
    paste0("--SAIGEOutputFile=", out),
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    quant_vr_model, quant_sparse_grm,
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01")
}
plink_geno <- c("--bedFile=./input/genotype_100markers.bed",
                "--bimFile=./input/genotype_100markers.bim",
                "--famFile=./input/genotype_100markers.fam",
                "--AlleleOrder=alt-first")
pgen_geno <- c("--pgenFile=./input/genotype_100markers.pgen",
               "--pvarFile=./input/genotype_100markers.pvar",
               "--psamFile=./input/genotype_100markers.psam")

test_that("P38: quantitative SKAT-O WITH a sparse GRM reproduces its golden", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- tempfile(fileext = ".txt")
  expect_cli_ok(run_saige_cli("step2_SPAtests.R",
                              c(quant_sparse_group_args(out, plink_geno), "--is_fastTest=FALSE")))
  gold <- golden_step2("plink_group_nofast_quant_sparse.txt")
  expect_output_matches_golden(gold, out)
  expect_output_matches_golden(paste0(gold, ".singleAssoc.txt"),
                               paste0(out, ".singleAssoc.txt"))
})

test_that("P39: quantitative BURDEN WITH a sparse GRM reproduces its golden", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- tempfile(fileext = ".txt")
  expect_cli_ok(run_saige_cli("step2_SPAtests.R",
                              c(quant_sparse_group_args(out, plink_geno),
                                "--r.corr=1", "--is_fastTest=FALSE")))
  expect_output_matches_golden(golden_step2("plink_group_burden_quant_sparse.txt"), out)
})

test_that("P40: quantitative single-variant WITH a sparse GRM reproduces its golden", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- tempfile(fileext = ".txt")
  expect_cli_ok(run_saige_cli("step2_SPAtests.R",
                              c(plink_geno, paste0("--SAIGEOutputFile=", out),
                                "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                                quant_vr_model, quant_sparse_grm, "--is_fastTest=FALSE")))
  expect_output_matches_golden(golden_step2("plink_marker_sparse_quant.txt"), out)
})

test_that("P41: the PGEN reader gives the same quantitative sparse-GRM result", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- tempfile(fileext = ".txt")
  expect_cli_ok(run_saige_cli("step2_SPAtests.R",
                              c(quant_sparse_group_args(out, pgen_geno), "--is_fastTest=FALSE")))
  gold <- golden_step2("pgen_group_nofast_quant_sparse.txt")
  expect_output_matches_golden(gold, out)
  expect_output_matches_golden(paste0(gold, ".singleAssoc.txt"),
                               paste0(out, ".singleAssoc.txt"))
})

test_that("the sparse-GRM runs are not copies of their no-GRM twins", {
  # The trap this guards: a new golden that reproduces an existing one gates
  # nothing. Measured when they were added -- the sparse GRM moves the
  # quantitative region p-values by up to 23% and BETA_Burden by up to 29%.
  skip_unless_level("standard")
  skip_if_no_goldens()
  for (pair in list(c("plink_group_nofast_quant.txt", "plink_group_nofast_quant_sparse.txt"),
                    c("plink_group_burden_quant.txt", "plink_group_burden_quant_sparse.txt"))) {
    a <- read_saige_table(golden_step2(pair[1]))
    b <- read_saige_table(golden_step2(pair[2]))
    expect_identical(nrow(a), nrow(b))
    x <- suppressWarnings(as.numeric(a$Pvalue_Burden))
    y <- suppressWarnings(as.numeric(b$Pvalue_Burden))
    ok <- is.finite(x) & is.finite(y) & x != 0
    expect_gt(max(abs(y[ok] / x[ok] - 1)), 1e-3)
  }
  # and the quantitative variance-ratio file really does carry `sparse` rows,
  # which is the whole reason these runs are possible
  vr <- read.table(paste0(quant_vr_step1, ".varianceRatio.txt"))
  expect_true(any(vr$V2 == "sparse"))
})

# --- LEDGER #86, characterization: asserts what the code does TODAY ----------
# When #86 is fixed this test changes with the ledger row: the saturated cells
# become p < 1, the component-test bound holds, and the Cauchy row goes back to
# being an ACAT.
test_that("LEDGER #86: SKAT-O p saturates at 1 and reroutes the Cauchy row to Bonferroni", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  d <- read_saige_table(golden_step2("plink_group_nofast_quant.txt"))
  cells <- d[d$Group != "Cauchy", ]
  p <- as.numeric(cells$Pvalue)
  pb <- as.numeric(cells$Pvalue_Burden)
  ps <- as.numeric(cells$Pvalue_SKAT)

  # (a) SKAT-O saturates: some cells report EXACTLY 1 ...
  expect_gt(sum(p == 1), 0)
  # ... and every saturated one has both component tests strictly below 1, so
  # the 1 is not merely inherited from a degenerate component.
  sat <- p == 1
  expect_true(all(pb[sat] < 1 & ps[sat] < 1))
  # ... and among the cells that actually RUN the SKAT-O combination they are
  # exactly the ones whose smaller component p is large. The all-ultra-rare
  # cells (Number_rare == 0) are excluded because they degenerate to a single
  # test -- Pvalue == Pvalue_Burden == Pvalue_SKAT, no rho grid, so nothing
  # there can saturate; asserting that separately is the control that the
  # ordering below is about the combination and not about p being large.
  comb <- as.integer(cells$Number_rare) > 0
  expect_true(all(!sat[!comb]))
  expect_gt(min(pmin(pb, ps)[sat & comb]), max(pmin(pb, ps)[!sat & comb]))
  # The threshold measured directly on the kernel is ~0.777; the ordering is
  # asserted instead of the constant, which belongs to the SKAT dependency.

  # (b) the invariant the suite already asserts on the binary golden
  # (test-invariants-goldens.R: SKAT-O lies between its component tests) is
  # VIOLATED here, on exactly the saturated rows.
  above <- p > pmax(pb, ps) * (1 + 1e-12)
  expect_identical(which(above), which(sat))

  # (c) a saturated cell reroutes its region's Cauchy row: the printed value is
  # min(1, n * min(p)) over the cells, NOT the ACAT of them.
  acat <- function(x) 0.5 - atan(mean(tan((0.5 - x) * pi))) / pi
  reg <- "GENE2"                       # the region where the two differ
  rc <- as.numeric(cells$Pvalue[cells$Region == reg])
  cauchy <- as.numeric(d$Pvalue[d$Region == reg & d$Group == "Cauchy"])
  expect_gt(sum(rc == 1), 0)
  expect_equal(cauchy, min(1, length(rc) * min(rc)), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(cauchy, acat(rc), tolerance = 1e-8)))
})
