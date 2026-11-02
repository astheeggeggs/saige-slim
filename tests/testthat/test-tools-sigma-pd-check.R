# harness/tools/aou_sigma_pd_check.R: a component above --maxDenseComponent is not
# evaluated, and a bound or a per-model verdict stated from the evaluated blocks alone
# describes part of the matrix (LEDGER #184, found on All of Us amr 2026-11-02: the
# 10,003-sample component held the -19.9 eigenvalue while the small blocks printed a
# bound of 1e15 and "within the tau bound"). Synthetic thinned GRMs, n = 400, diagonal 1:
#   chain    samples 1..300 at 0.25 between neighbours: positive definite (eigenvalues
#            1 + 0.5 cos(k pi / 301), smallest 0.5), the large component
#   small    a star of 4 at 0.9 (hub 310, leaves 311..313): smallest eigenvalue
#            1 - 0.9 sqrt(3) = -0.5588, a slightly indefinite SMALL block (bound 1.789)
#   planted  a star of 4 at 100 INSIDE the chain (hub 150, leaves 301..303): the chain's
#            component is now 303 samples and indefinite (smallest eigenvalue -172.2, bound
#            0.0058), by interlacing whatever the chain adds
#   families (320,321) 0.5; (330,331) 0.5, (330,332) 0.25, (331,332) 0.25: positive definite
# Fixtures: clean = chain + families; small = clean + the small star; planted = small + the
# star in the chain. A quantitative model at tau2/tau1 = 0.1 sits inside the small star's
# bound and outside the planted one's. The reference for every bound is a dense eigen of
# the whole 400 x 400 computed here, independently of the script's block decomposition.
# Drilled 2026-11-02: the unpatched script on `planted` at cap 200 printed the bound 1.789
# and judged the model positive definite (margin 0.944); the fix refuses both.

# SIGMA_PD_SCRIPT points the file at another copy of the script: the drill against the
# unpatched one (git show <rev>:harness/tools/aou_sigma_pd_check.R > old.R)
sigma_script <- Sys.getenv("SIGMA_PD_SCRIPT",
                           unset = file.path(saige_repo, "harness", "tools",
                                             "aou_sigma_pd_check.R"))

sigma_fixture <- function(dir, kind) {
  n <- 400L
  i <- 1:299
  j <- 2:300
  x <- rep(0.25, 299)
  if (kind == "planted") {
    i <- c(i, rep(150L, 3))
    j <- c(j, 301:303)
    x <- c(x, rep(100, 3))
  }
  if (kind != "clean") {
    i <- c(i, rep(310L, 3))
    j <- c(j, 311:313)
    x <- c(x, rep(0.9, 3))
  }
  i <- c(i, 320L, 330L, 330L, 331L)
  j <- c(j, 321L, 331L, 332L, 332L)
  x <- c(x, 0.5, 0.5, 0.25, 0.25)
  K <- Matrix::sparseMatrix(i = c(i, 1:n), j = c(j, 1:n), x = c(x, rep(1, n)),
                            dims = c(n, n), symmetric = TRUE)
  f <- file.path(dir, paste0(kind, ".mtx"))
  Matrix::writeMM(K, f)
  writeLines(sprintf("S%04d", 1:n), paste0(f, ".sampleIDs.txt"))
  f
}

sigma_reference_bound <- function(mtx) {
  ev <- eigen(as.matrix(Matrix::readMM(mtx)), symmetric = TRUE, only.values = TRUE)$values
  c(min = min(ev), bound = 1 / abs(min(ev)))
}

sigma_models <- function(dir) {
  d <- file.path(dir, "models")
  dir.create(d, showWarnings = FALSE)
  modglmm <- list(theta = c(1, 0.1), traitType = "quantitative", sampleID = sprintf("S%04d", 1:400))
  save(modglmm, file = file.path(d, "Q_ratio0.1.rda"))
  d
}

sigma_check <- function(mtx, cap, model_dir = NULL, label = "t") {
  rep <- tempfile(fileext = ".tsv")
  log <- tempfile(fileext = ".log")
  args <- c(sigma_script, "--grmOnly", paste0("--sparseGRMFile=", mtx),
            paste0("--sparseGRMSampleIDFile=", mtx, ".sampleIDs.txt"),
            "--relatednessCutoff=0.05", paste0("--maxDenseComponent=", cap),
            paste0("--reportFile=", rep), paste0("--label=", label))
  if (!is.null(model_dir)) args <- c(args, paste0("--modelDir=", model_dir))
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), args,
                                 stdout = log, stderr = log))
  testthat::expect_identical(as.integer(st), 0L,
                             info = paste(readLines(log, warn = FALSE), collapse = "\n"))
  tab <- utils::read.table(rep, header = TRUE, sep = "\t", colClasses = "character",
                           stringsAsFactors = FALSE)
  structure(stats::setNames(tab$value, tab$key), log = readLines(log, warn = FALSE))
}

sigma_dir <- local({
  d <- tempfile("sigma-pd-")
  dir.create(d)
  d
})
fx <- vapply(c("clean", "small", "planted"), function(k) sigma_fixture(sigma_dir, k), character(1))
ref <- lapply(fx, sigma_reference_bound)
md <- sigma_models(sigma_dir)

test_that("the fixtures are what the header says (the reference, not the script)", {
  expect_equal(unname(ref$clean["min"]), 0.5, tolerance = 1e-6)
  expect_equal(unname(ref$small["min"]), 1 - 0.9 * sqrt(3), tolerance = 1e-6)
  expect_lt(unname(ref$planted["min"]), 1 - 100 * sqrt(3) + 1e-6)   # interlacing: the star's
  expect_lt(unname(ref$planted["bound"]), 0.01)
})

test_that("a planted indefinite block ABOVE the cap: no bound, no verdict, no model judged", {
  r <- sigma_check(fx[["planted"]], cap = 200, model_dir = md)
  expect_identical(r[["components_not_evaluated"]], "1")
  expect_identical(r[["samples_not_evaluated"]], "303")
  expect_false("tau2_over_tau1_bound_quantitative" %in% names(r))
  expect_false("tau2_bound_binary" %in% names(r))
  expect_false(any(startsWith(names(r), "model_")))
  expect_match(r[["verdict"]],
               "^GRM indefinite: bound unknown, 1 component\\(s\\) \\(303 samples\\) above")
  expect_match(r[["models_checked"]], "^0 ")
  expect_true(any(grepl("above --maxDenseComponent 200: NOT evaluated", attr(r, "log"),
                        fixed = TRUE)))
  # the same without --modelDir
  r2 <- sigma_check(fx[["planted"]], cap = 200)
  expect_match(r2[["verdict"]], "^GRM indefinite: bound unknown")
  expect_false("tau2_over_tau1_bound_quantitative" %in% names(r2))
})

test_that("the same block once the cap covers it: the reference bound, and the model refused", {
  r <- sigma_check(fx[["planted"]], cap = 400, model_dir = md)
  expect_identical(r[["components_not_evaluated"]], "0")
  expect_equal(as.numeric(r[["min_eigenvalue_sigma"]]), unname(ref$planted["min"]),
               tolerance = 1e-5)
  expect_equal(as.numeric(r[["tau2_over_tau1_bound_quantitative"]]), unname(ref$planted["bound"]),
               tolerance = 1e-5)
  expect_match(r[["model_Q_ratio0.1"]], "NOT positive definite")
  expect_identical(r[["verdict"]], "Sigma is NOT positive definite for some model in --modelDir")
  expect_false(any(grepl("^dense block of ", attr(r, "log"))))   # 303 samples: below the log line
})

test_that("a small indefinite block with a large PD component skipped is ALSO undecided", {
  # the AoU shape in miniature: the unpatched script printed the small block's bound here
  r <- sigma_check(fx[["small"]], cap = 200, model_dir = md)
  expect_identical(r[["components_not_evaluated"]], "1")
  expect_identical(r[["samples_not_evaluated"]], "300")
  expect_false("tau2_over_tau1_bound_quantitative" %in% names(r))
  expect_match(r[["verdict"]], "^GRM indefinite: bound unknown")
  # and decided once evaluated: the small star's bound, the model inside it
  r <- sigma_check(fx[["small"]], cap = 400, model_dir = md)
  expect_identical(r[["components_not_evaluated"]], "0")
  expect_equal(as.numeric(r[["tau2_over_tau1_bound_quantitative"]]), unname(ref$small["bound"]),
               tolerance = 1e-5)
  expect_match(r[["model_Q_ratio0.1"]], ": positive definite$")
  expect_identical(r[["verdict"]],
                   "GRM indefinite: Sigma positive definite for every model in --modelDir")
})

test_that("an unplanted large PSD component still reads positive definite, cap above or below it", {
  for (cap in c(200, 400)) {
    r <- sigma_check(fx[["clean"]], cap = cap)
    expect_identical(r[["sparse_cholesky_succeeds"]], "TRUE")
    expect_identical(r[["verdict"]], "positive definite")
    expect_identical(r[["components_not_evaluated"]], if (cap == 200) "1" else "0")
  }
})

test_that("the bundled sparse GRM and binary model read as before, with the three new rows", {
  g <- file.path(saige_extdata, "output",
                 "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
  skip_if_not(file.exists(g))
  rep <- tempfile(fileext = ".tsv")
  model <- file.path(saige_extdata, "output", "example_binary_sparseGRM.rda")
  args <- c(sigma_script, paste0("--modelFile=", model), paste0("--sparseGRMFile=", g),
            paste0("--sparseGRMSampleIDFile=", g, ".sampleIDs.txt"), "--relatednessCutoff=0.125",
            paste0("--reportFile=", rep), "--label=b")
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), args,
                                 stdout = FALSE, stderr = FALSE))
  expect_identical(as.integer(st), 0L)
  tab <- utils::read.table(rep, header = TRUE, sep = "\t", colClasses = "character")
  r <- stats::setNames(tab$value, tab$key)
  expect_identical(r[["verdict"]], "positive definite")
  expect_identical(r[["components_not_evaluated"]], "0")
  expect_identical(r[["samples_not_evaluated"]], "0")
  expect_true("largest_dense_block_gb" %in% names(r))
})
