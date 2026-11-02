# Step 0's pair-major refine (refineKinPairMajor in src/SAIGE_fitGLMM_fast.cpp,
# LEDGER #185, 2026-11-02) against the marker-major loop it replaces
# (SAIGE_REFINE_PAIRMAJOR=0).
#
# The claim under test is BIT IDENTITY, by construction: per screened pair the
# refined kinship is the same sequence of float additions (kin = table + kin,
# marker by marker, the blocks taken in order), reading the per-marker 3 x 3
# table setSparseKinLookUpArr fills and the genotypes from the same packed
# store. What changes is only the loop order (pairs outer, markers inner, a
# block of markers transposed once into sample-major bytes) -- the loop that
# was 3,121 s of a 4,199 s All of Us eur step 0 (33.8M screened pairs x 99,900
# markers, 64 cores).
#
# Controls:
#   - every C(N,2) pair at cutoff -1e30: the two loops give identical values,
#     for block lengths that do and do not divide M, block 1, block > M, and
#     1 and 3 threads; the run REPORTS which loop ran (a knob that fails to
#     arrive would otherwise pass as green, the LEDGER #30 shape);
#   - the planted parent-offspring and duplicate pairs refine to ~0.5 and ~1
#     (positive control that the values are kinships, not zeros);
#   - the CLI on the bundled fixture reproduces the committed step-0 golden,
#     through both readers, at a block length that does not divide M.

test_that("pair-major refine is bit-identical to the marker-major loop, any block, any threads", {
  skip_unless_level("standard")
  dir <- tempfile()
  dir.create(dir)
  prefix <- file.path(dir, "ref")
  fx <- make_related_geno(400L, 150L, 6L, 4L, 11L)
  write_plink_fixture_b(prefix, fx$g)
  n <- ncol(fx$g)
  old <- run_screen_b(prefix, n, "pairwise", cutoff = -1e30, refine = "markermajor")
  expect_identical(length(old$k), as.integer(choose(n, 2)))
  expect_true(any(grepl("^refine: marker-major", old$log)))
  expect_false(any(grepl("^refine: pair-major", old$log)))
  # positive control: planted pairs refine to kinships
  pl <- match(fx$planted, old$keys)
  expect_false(anyNA(pl))
  expect_true(all(old$k[pl] > 0.35))
  for (b in list(NULL, 1L, 7L, 333L, 100000L)) {
    for (t in c(1L, 3L)) {
      new <- run_screen_b(prefix, n, "pairwise", cutoff = -1e30, refine = "pairmajor",
                          refine_block = b, threads = t)
      expect_identical(new$keys, old$keys)
      expect_identical(new$k, old$k)
      expect_true(any(grepl("^refine: pair-major, [0-9]+ markers in [0-9]+ block", new$log)))
    }
  }
})

test_that("pair-major refine reproduces the step-0 golden byte for byte, both readers", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  suffix <- "_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
  withr::local_envvar(SAIGE_REFINE_BLOCK = "333")
  for (reader in c("plinkFile", "plink2File")) {
    out <- file.path(tempdir(), paste0("s0_pairmajor_", reader))
    st <- run_saige_cli("createSparseGRM.R", c(
      paste0("--", reader, "=input/plinkforGRM_1000samples_10kMarkers"),
      "--nThreads=2", paste0("--outputPrefix=", out),
      "--numRandomMarkerforSparseKin=5000", "--relatednessCutoff=0.05",
      "--screenMethod=blocked"
    ))
    expect_cli_ok(st)
    golden <- file.path(saige_goldens, "step0",
                        if (reader == "plinkFile") "test_sparseGRM" else "test_sparseGRM_pgen")
    expect_files_byte_identical(paste0(golden, suffix), paste0(out, suffix))
    expect_files_byte_identical(paste0(golden, suffix, ".sampleIDs.txt"),
                                paste0(out, suffix, ".sampleIDs.txt"))
  }
})
