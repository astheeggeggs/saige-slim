# --is_output_nEff (2026-10-16): Kish's effective number of carriers of each
# cell's burden vector, (sum x)^2 / sum x^2 over the per-sample weighted sums
# the burden test saw (the ultra-rare pseudo-marker included), as the column
# nEff_Burden after SE_Burden (and after the log10 columns when both are on),
# on both region paths, NA on the pooled row. Why: a steeply weighted set of
# many singletons can rest on two or three people, and on a quantitative trait
# that is an outlier person rather than a gene (AoU, runbook 19.8: three eur
# Height FlexRV hits on fewer than five effective carriers). The column lets a
# meta-analysis filter such cells. It changes no p-value.
#   (1) off: a run with the flag is the run without it in every pre-existing
#       column, both paths, and the run REPORTS the flag;
#   (2) on: the column sits where the schema says, the pooled row is NA, and
#       every cell obeys the bounds a Kish count must: 0 < nEff <= MAC (the
#       count cannot exceed the carriers, who cannot exceed the alleles), and
#       on an UNWEIGHTED set every carrier's value is 1 or 2, so nEff >= MAC/2;
#   (3) the C++ arithmetic against closed forms (equal entries give n, one
#       dominant entry gives ~1, an all-zero vector gives 0);
#   (4) the two goldens P81/P82 (test-step2-goldens.R) pin the bytes.

neff_run <- function(path, flag, out = tempfile(fileext = ".txt"), extra = character()) {
  m <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
  sg_name <- "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
  sg <- file.path(saige_extdata, "output", sg_name)
  common <- c("--bedFile=./input/genotype_100markers.bed",
              "--bimFile=./input/genotype_100markers.bim",
              "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
              "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
              paste0("--GMMATmodelFile=", m, ".rda"),
              paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
              paste0("--sparseGRMFile=", sg),
              paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
              paste0("--is_output_nEff=", flag),
              paste0("--SAIGEOutputFile=", out), extra)
  args <- if (path == "skato") {
    c(common, "--groupFile=./input/group_new_snpid.txt",
      "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
      "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=TRUE")
  } else {
    c(common, "--groupFile=./input/group_new_snpid_flexrv.txt",
      "--flexRV_score=REVEL", "--flexRV_maxMAF=0.01",
      "--annotation_in_groupTest=lof:missense", "--r.corr=1", "--is_fastTest=FALSE")
  }
  st <- run_saige_cli("step2_SPAtests.R", args)
  expect_cli_ok(st)
  list(out = out, log = readLines(attr(st, "log"), warn = FALSE))
}

read_rows <- function(f) read.delim(f, check.names = FALSE, colClasses = "character")

test_that("the flag is additive on both paths, and the run reports it", {
  skip_unless_level("standard")
  for (path in c("skato", "burden")) {
    off <- neff_run(path, "FALSE")
    on <- neff_run(path, "TRUE")
    a <- read_rows(on$out)
    b <- read_rows(off$out)
    expect_identical(names(a)[names(a) != "nEff_Burden"], names(b), info = path)
    expect_identical(a[, names(b)], b, info = path)
    expect_identical(names(a)[match("SE_Burden", names(a)) + 1], "nEff_Burden", info = path)
    expect_true(any(grepl("nEff_Burden column (--is_output_nEff) is on", on$log, fixed = TRUE)),
                info = path)
    expect_false(any(grepl("--is_output_nEff) is on", off$log, fixed = TRUE)), info = path)
  }
})

test_that("the column is a Kish count: NA on the pooled row, 0 < nEff <= MAC on every cell", {
  skip_unless_level("standard")
  for (path in c("skato", "burden")) {
    d <- read_rows(neff_run(path, "TRUE")$out)
    pooled <- d$Group == "Cauchy"
    expect_true(all(is.na(d$nEff_Burden[pooled])), info = path) # written "NA", read back as NA
    ne <- as.numeric(d$nEff_Burden[!pooled])
    mac <- as.numeric(d$MAC[!pooled])
    expect_true(all(is.finite(ne)), info = path)
    expect_true(all(ne[mac > 0] > 0), info = path)
    expect_true(all(ne <= mac + 1e-9), info = path)
  }
})

test_that("on an unweighted set every carrier holds 1 or 2, so MAC/2 <= nEff <= MAC", {
  skip_unless_level("standard")
  # --weights.beta=1,1 makes dbeta(MAF, 1, 1) = 1 for every marker: the burden
  # vector is the plain per-sample sum of genotypes over the cell's markers
  d <- read_rows(neff_run("skato", "TRUE", extra = "--weights.beta=1,1")$out)
  cells <- d$Group != "Cauchy"
  ne <- as.numeric(d$nEff_Burden[cells])
  mac <- as.numeric(d$MAC[cells])
  expect_true(all(ne >= mac / 2 - 1e-9 & ne <= mac + 1e-9))
  # and a cell whose carriers are all single hets of one allele each is exact: nEff == MAC
  n_var <- as.numeric(d$Number_rare) + as.numeric(d$Number_ultra_rare)
  singles <- cells & as.numeric(d$MAC) == n_var
  if (any(singles)) {
    expect_equal(as.numeric(d$nEff_Burden[singles]), as.numeric(d$MAC[singles]), tolerance = 1e-9)
  }
})

test_that("the C++ Kish arithmetic matches its closed forms", {
  expect_equal(kishNEffInCPP(rep(1, 7)), 7)
  expect_equal(kishNEffInCPP(rep(0.3, 12)), 12)
  expect_equal(kishNEffInCPP(c(1, 1e-6, 1e-6)), (1 + 2e-6)^2 / (1 + 2e-12), tolerance = 1e-12)
  expect_lt(kishNEffInCPP(c(1, 1e-6, 1e-6)), 1.00001)
  expect_identical(kishNEffInCPP(numeric(0)), 0)
  expect_identical(kishNEffInCPP(c(0, 0, 0)), 0)
  x <- c(0.2, 0.5, 1, 2)
  expect_equal(kishNEffInCPP(x), sum(x)^2 / sum(x^2), tolerance = 1e-12)
})
