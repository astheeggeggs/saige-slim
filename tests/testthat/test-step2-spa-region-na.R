# LEDGER #181 (2026-10-16): the saddlepoint the region path runs from R
# (get_newPhi_scaleFactor_traitType -> SPA_pval) handed a NaN into
# SPA_ER_kernel_related_Phiadj_fast_new's `if (VarQ_2 == 0)`, which errored
# OUTSIDE the per-cell catch and killed the chromosome's R process (afr HTN,
# four production jobs, delivered as SUCCESS with chrStatus the only record).
# Now a non-finite rescale p is counted, named, dumped under SAIGE_SPA_DUMP_DIR
# and replaced by the normal-approximation p with the convergence flag FALSE.
# No bundled cell produces the NaN, so the net is exercised by a plant,
# SAIGE_SPA_PLANT_NA_R=1 (every rescale saddlepoint returns NA):
#   negative control  P2 as shipped: the golden, the line reads 0, no sidecar;
#   positive control  P2 planted: the run COMPLETES, the line counts the
#                     rescales that reached the saddlepoint, the sidecar names
#                     them with the normal p that was reported, every dump has
#                     the vector, and every cell's p-values stay finite.
skip_if_no_goldens()

na_args <- function(out) {
  m <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
  sg_name <- "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  sg <- file.path(saige_extdata, "output", sg_name)
  c("--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
    paste0("--SAIGEOutputFile=", out),
    "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
    paste0("--GMMATmodelFile=", m, ".rda"),
    paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
    paste0("--sparseGRMFile=", sg),
    paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
    "--groupFile=./input/group_new_snpid.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=TRUE")
}
golden_p2 <- file.path(saige_goldens, "step2", "plink_group_fastTest.txt")
na_line <- function(log) grep("^Region saddlepoint NA: ", log, value = TRUE)

test_that("as shipped: the golden, the line reads 0 and there is no sidecar", {
  skip_unless_level("standard")
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", na_args(out))
  expect_cli_ok(st)
  expect_output_matches_golden(golden_p2, out)
  log <- readLines(attr(st, "log"), warn = FALSE)
  # on the bundled model every rescale saddlepoint converges: 0, and no sidecar
  expect_match(na_line(log), "^Region saddlepoint NA: 0 burden rescale saddlepoints failed")
  expect_false(file.exists(paste0(out, ".spaRegionNA.txt")))
})

test_that("planted: the run completes, every rescale is counted, named and dumped, no p is NA", {
  skip_unless_level("standard")
  out <- tempfile(fileext = ".txt")
  dir <- tempfile("spadump_")
  dir.create(dir)
  st <- run_saige_cli("step2_SPAtests.R", na_args(out),
                      env = c(SAIGE_SPA_PLANT_NA_R = "1", SAIGE_SPA_DUMP_DIR = dir))
  expect_cli_ok(st)
  log <- readLines(attr(st, "log"), warn = FALSE)
  line <- na_line(log)
  n <- as.integer(sub("^Region saddlepoint NA: ([0-9]+) burden rescale.*$", "\\1", line))
  expect_match(line, paste0(" ", n, " of them returned a non-finite value"))
  expect_gt(n, 0)
  side <- paste0(out, ".spaRegionNA.txt")
  expect_true(file.exists(side))
  tab <- read.delim(side, stringsAsFactors = FALSE)
  expect_identical(nrow(tab), n)
  expect_true(all(tab$Kind == "region_rescale"))
  expect_true(all(tab$Reason == "nonfinite_reached_R"))
  expect_true(all(tab$Region %in% c("GENE1", "GENE2")))
  # the plant nulls the value AFTER the kernel, so the recorded kernel value is the converged one
  expect_true(all(is.finite(tab$kernel_value)))
  expect_true(all(is.finite(tab$log_p_normal_reported)))
  # every named cell has its dump, and the dump holds the model's vector
  expect_true(all(nzchar(tab$dump)) && all(file.exists(tab$dump)))
  d <- read.delim(tab$dump[1], comment.char = "#")
  expect_identical(names(d), c("mu", "g"))
  expect_identical(nrow(d), 1000L)
  hdr <- readLines(tab$dump[1], n = 1)
  expect_match(hdr, "^# kernel=SPA_pval\\(region rescale\\) q=.* n=1000$")
  # the results: same cell grid as the golden, every p-value finite (the
  # normal p was used for the rescale; the SKAT-O kernel still ran)
  gold <- read_saige_table(golden_p2)
  new <- read_saige_table(out)
  expect_identical(dim(new), dim(gold))
  expect_identical(new[, c("Region", "Group", "max_MAF")], gold[, c("Region", "Group", "max_MAF")])
  pv <- suppressWarnings(as.numeric(new$Pvalue[new$Group != "Cauchy"]))
  expect_true(all(is.finite(pv)))
  # the plant moves only cells whose rescale reached the saddlepoint (|z| >= 2)
  moved <- which(new$Pvalue != gold$Pvalue)
  expect_lte(length(moved), nrow(new))
})
