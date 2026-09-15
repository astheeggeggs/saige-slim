# LEDGER #153 (FIXED 2026-09-15): getNeff solves Nglmm on the sparse GRM and
# refuses to run without one.
#
# Nglmm is 1' K^-1 1, the effective number of unrelated samples. Before the
# fix a sparse GRM given WITHOUT --useSparseGRMtoFitNULL=TRUE was never loaded
# into C++, so the PCG ran on the dense 128,868-marker GRM with no residual
# term, hit maxiterPCG, and its partial sum was printed as Nglmm: 9.672876e+15
# for N = 1,000. The same run with the flag printed 402.4175, which is what
# solve(K, 1) summed in R gives from the same .mtx (LEDGER row). Now the file
# alone is enough, the flag is implied, and no file is a stop.
#
# The number is checked two ways: the run's printed Nglmm against the dense
# closed form computed here from the .mtx, at 1e-6 relative (the print has
# seven significant digits). Tier: standard (one extractNglmm run + a stop).

l153_cli <- function(...) get("run_saige_cli")(...)

l153_sg <- "output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"

l153_args <- function(with_grm) {
  c(
    "--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr",
    if (with_grm) c(paste0("--sparseGRMFile=", l153_sg),
                    paste0("--sparseGRMSampleIDFile=", l153_sg, ".sampleIDs.txt")),
    "--phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt",
    "--phenoCol=y_quantitative", "--covarColList=x1,x2", "--sampleIDColinphenoFile=IID",
    "--traitType=quantitative", paste0("--outputPrefix=", tempfile("nglmm153_"))
  )
}

test_that("getNeff with a sparse GRM and no flag solves on the sparse GRM (LEDGER #153)", {
  skip_unless_level("standard")
  st <- l153_cli("extractNglmm.R", l153_args(TRUE))
  lg <- readLines(attr(st, "log"), warn = FALSE)
  expect_identical(as.integer(st), 0L)
  expect_true(any(grepl("--useSparseGRMtoFitNULL=TRUE is implied", lg, fixed = TRUE)))
  expect_true(any(grepl("use sparse kinship to fit the model", lg, fixed = TRUE)))
  expect_false(any(grepl("pcg did not converge", lg, fixed = TRUE)))
  ng <- grep("^Nglmm ", lg, value = TRUE)
  expect_length(ng, 1L)
  nglmm <- as.numeric(sub("^Nglmm ", "", ng))
  # the closed form from the same matrix, the second way
  K <- Matrix::forceSymmetric(
    methods::as(Matrix::readMM(file.path(saige_extdata, l153_sg)), "CsparseMatrix"), uplo = "L")
  closed <- sum(solve(as.matrix(K), rep(1, nrow(K))))
  expect_gt(closed, 1); expect_lt(closed, nrow(K))
  expect_equal(nglmm, closed, tolerance = 1e-6)
  expect_lt(nglmm, 1e6) # the number the row opened on was 9.67e15
})

test_that("getNeff without a sparse GRM stops naming the row (LEDGER #153)", {
  skip_unless_level("standard")
  st <- l153_cli("extractNglmm.R", l153_args(FALSE))
  lg <- readLines(attr(st, "log"), warn = FALSE)
  expect_false(as.integer(st) == 0L)
  expect_true(any(grepl("getNeff needs --sparseGRMFile", lg, fixed = TRUE)))
  expect_false(any(grepl("^Nglmm ", lg)))
})
