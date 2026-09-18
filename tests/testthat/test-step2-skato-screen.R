# LEDGER #158: the SKAT-O outer-integral screen (R/SkatoScreen.R).
#
# Three things are pinned. (1) With the screen OFF (--SKATO_screen_threshold=0)
# the run is the upstream kernel byte for byte: golden P2 reproduced, no
# sidecar, and the log SAYS the screen is off. (2) With the screen ON, only
# `Pvalue` may differ from the off run, only on cells the sidecar marks
# "liu-screen", and by under 15%; Pvalue_SKAT, Pvalue_Burden, BETA_Burden,
# SE_Burden and every other column are byte-identical, the miss count is zero
# (no screened cell with a Davies p below threshold/10), and the log reports
# the threshold and the counts, which add up to the sidecar. (3) The kernel
# reimplementation reproduces SKAT's own answer exactly on random cells when
# forced down the Davies path, and its per-rho layer equals SKAT's on every
# cell -- the control on the transcription.
#
# The on/off comparison is against the OFF run of the same session, not the
# committed golden, so it survives the goldens being re-baselined to the
# screened values.

screen_run <- function(threshold, out = tempfile(fileext = ".txt")) {
  m <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
  sg <- file.path(saige_extdata, "output",
    "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
  args <- c("--bedFile=./input/genotype_100markers.bed",
            "--bimFile=./input/genotype_100markers.bim",
            "--famFile=./input/genotype_100markers.fam", "--AlleleOrder=alt-first",
            "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
            paste0("--GMMATmodelFile=", m, ".rda"),
            paste0("--varianceRatioFile=", m, ".varianceRatio.txt"),
            paste0("--sparseGRMFile=", sg),
            paste0("--sparseGRMSampleIDFile=", sg, ".sampleIDs.txt"),
            "--groupFile=./input/group_new_snpid.txt",
            "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
            "--maxMAF_in_groupTest=0.0001,0.001,0.01", "--is_fastTest=TRUE",
            paste0("--SKATO_screen_threshold=", threshold),
            paste0("--SAIGEOutputFile=", out))
  st <- run_saige_cli("step2_SPAtests.R", args)
  list(status = as.integer(st), log = readLines(attr(st, "log"), warn = FALSE), out = out)
}

read_region <- function(path) {
  utils::read.table(path, header = TRUE, sep = "\t", colClasses = "character",
                    quote = "", comment.char = "", na.strings = character(0))
}

test_that("screen OFF is the upstream kernel byte for byte, and says so", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- screen_run(0)
  expect_cli_ok(r$status)
  expect_true(any(grepl("SKAT-O screen: off", r$log, fixed = TRUE)))
  expect_false(file.exists(paste0(r$out, ".skatoMethod.txt")))
  expect_files_byte_identical(file.path(saige_goldens, "step2", "plink_group_fastTest_screenoff.txt"), r$out)
})

test_that("screen ON touches only Pvalue, only on screened cells, by under 15%, with zero misses", {
  skip_unless_level("standard")
  off <- screen_run(0)
  on <- screen_run(0.05)
  expect_cli_ok(off$status)
  expect_cli_ok(on$status)
  line <- grep("^SKAT-O screen: threshold 0.05", on$log, value = TRUE)
  expect_length(line, 1)
  side <- paste0(on$out, ".skatoMethod.txt")
  expect_true(file.exists(side))
  sc <- read_region(side)
  nScreened <- as.integer(sub(".*; ([0-9]+) cell\\(s\\) took the Liu-form.*", "\\1", line))
  nDavies <- as.integer(sub(".* ([0-9]+) the Davies integral.*", "\\1", line))
  expect_equal(nScreened, sum(sc$Method == "liu-screen"))
  expect_equal(nDavies, sum(sc$Method == "davies"))
  expect_gt(nScreened, 0)
  a <- read_region(off$out)
  b <- read_region(on$out)
  expect_identical(dim(a), dim(b))
  expect_identical(names(a), names(b))
  cells <- a$Group != "Cauchy"
  for (col in setdiff(names(a), "Pvalue")) {
    expect_identical(a[[col]][cells], b[[col]][cells], label = paste("column", col))
  }
  key <- function(d) paste(d$Region, d$Group, d$max_MAF)
  sc <- sc[sc$Test == "marginal", ]
  screened <- key(a) %in% key(sc)[sc$Method == "liu-screen"]
  expect_identical(a$Pvalue[cells & !screened], b$Pvalue[cells & !screened])
  pa <- as.numeric(a$Pvalue[cells & screened])
  pb <- as.numeric(b$Pvalue[cells & screened])
  expect_true(all(abs(pb / pa - 1) < 0.15))
  expect_true(all(pb > 0.05))
  # the miss count: no screened cell whose Davies p is below threshold / 10
  expect_equal(sum(pa < 0.005), 0)
  # the Cauchy rows are the only other rows that may move, through the cells
  expect_true(all(abs(as.numeric(b$Pvalue[!cells]) / as.numeric(a$Pvalue[!cells]) - 1) < 0.15))
})

test_that("transcribed kernel reproduces SKAT exactly (control); screens only the outer integral", {
  set.seed(158)
  rho <- c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1)
  for (q in c(2, 5, 12)) {
    G <- matrix(rnorm(40 * q), 40, q)
    Phi <- crossprod(G) / 40
    Score <- as.vector(t(chol(Phi)) %*% rnorm(q)) * 1.5
    ref <- SKAT:::Met_SKAT_Get_Pvalue(Score = Score, Phi = Phi, r.corr = rho,
                                      method = "optimal.adj", Score.Resampling = NULL)
    # forced down the Davies path: identical object
    skato_screen_env$threshold <- 1
    got <- skatoScreenKernel(Score, Phi, rho)
    expect_identical(attr(got, "skato_method"), "davies")
    attr(got, "skato_method") <- NULL
    attr(got, "skato_screen_p") <- NULL
    expect_identical(got, ref)
    # forced down the screen: per-rho layer identical, p.value the Liu form
    skato_screen_env$threshold <- .Machine$double.eps
    got2 <- skatoScreenKernel(Score, Phi, rho)
    expect_identical(attr(got2, "skato_method"), "liu-screen")
    expect_identical(got2$param$p.val.each, ref$param$p.val.each)
    expect_identical(got2$param$rho, ref$param$rho)
    expect_identical(got2$param$minp, ref$param$minp)
    expect_true(abs(got2$p.value / ref$p.value - 1) < 0.15)
  }
  skato_screen_env$threshold <- 0
})
