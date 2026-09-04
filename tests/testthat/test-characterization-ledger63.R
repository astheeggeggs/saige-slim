# Characterization test for LEDGER.md #63 (OPEN): the sparse-GRM relatedness
# screen enumerates C(N,2) - 1 sample pairs, not C(N,2).
#
# `findIndiceRelatedSample` hands a HALF-OPEN parallelFor the count
# `Ntotal*(Ntotal-1)/2 - 1`, and the k -> (i,j) triangle fold inside the worker
# is an exact bijection onto all C(N,2) pairs over k = 0..C(N,2)-1. So the `- 1`
# drops one REAL pair -- the image of the last k -- and `geno.indiceVec` is the
# only pair list `refineKin` ever re-scores, so a pair never screened can never
# reach the .mtx.
#
# This test asserts TODAY's behaviour. When #63 is fixed, both the ledger row
# and this file change together: the counts become choose(N, 2) and the
# `missing` sets become empty.
#
# The fixture is derived at run time, and it is tiny (N = 10 and 11) so that
# the screen can be run with a cutoff low enough to retain EVERY pair it visits
# -- which is what makes the returned pair list a direct readout of the loop
# bound rather than of the relatedness filter. No committed fixture can serve:
# step 0's is N = 1,000, where retaining all 499,500 pairs is not cheap and
# where the dropped pair happens to be unrelated (kinship -0.0087) and so is
# invisible in the output.

# Write a plink1 .bed/.bim/.fam triple. `geno` is M x N, entries 0/1/2 counting
# copies of A1 (bim column 5); variant-major, 2 bits per sample, LSB-first,
# codes 00 = 2 copies, 10 = 1, 11 = 0 (01 = missing, unused here).
write_plink_fixture <- function(prefix, geno) {
  m <- nrow(geno)
  n <- ncol(geno)
  writeLines(sprintf("1\tsnp%d\t0\t%d\tA\tT", seq_len(m), seq_len(m)),
             paste0(prefix, ".bim"))
  writeLines(sprintf("f%d\ts%d\t0\t0\t0\t-9", seq_len(n), seq_len(n)),
             paste0(prefix, ".fam"))
  code <- c(3L, 2L, 0L)                      # geno 0/1/2 -> 2-bit code
  nb <- ceiling(n / 4)
  con <- file(paste0(prefix, ".bed"), "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  for (i in seq_len(m)) {
    cd <- c(code[geno[i, ] + 1L], rep(0L, nb * 4L - n))
    bytes <- vapply(seq_len(nb), function(b) {
      u <- cd[(4L * (b - 1L) + 1L):(4L * b)]
      as.raw(u[1] + 4L * u[2] + 16L * u[3] + 64L * u[4])
    }, raw(1))
    writeBin(bytes, con)
  }
  invisible(prefix)
}

# Deterministic genotypes, every marker polymorphic so every marker passes the
# MAF/missingness QC and the screen sees all M of them.
make_geno <- function(m, n, seed) {
  set.seed(seed)
  g <- matrix(rbinom(m * n, 2, 0.3), nrow = m, ncol = n)
  for (i in seq_len(m)) {
    if (sum(g[i, ]) == 0) g[i, 1] <- 1L
    if (sum(g[i, ]) == 2 * n) g[i, 1] <- 1L
  }
  g
}

# Ingest the fixture and run the screen in a child process: setgeno and the
# screen both write C++ globals, so they cannot share this session.
run_screen <- function(prefix, n, cutoff = -1e6) {
  out <- tempfile(fileext = ".rds")
  log <- tempfile(fileext = ".log")
  script <- tempfile(fileext = ".R")
  writeLines(c(
    "suppressPackageStartupMessages(library(SAIGE))",
    sprintf('prefix <- "%s"; out <- "%s"; n <- %d; cutoff <- %g',
            prefix, out, n, cutoff),
    "setminMAFforGRM(0.01); setmaxMissingRateforGRM(0.15)",
    'setgeno(paste0(prefix, ".bed"), paste0(prefix, ".bim"),',
    '        paste0(prefix, ".fam"), seq_len(n), rep(TRUE, n), 2, FALSE)',
    # Every QC-passing marker is a sub-marker: no sample(), so no RNG.
    "mq <- length(getAlleleFreqVec())",
    "setRelatednessCutoff(cutoff)",
    "setSubMarkerIndex(seq_len(mq) - 1L)",
    "Get_MultiMarkersBySample_StdGeno_Mat()",
    "std <- vapply(seq_len(n) - 1L, getColfromStdGenoMultiMarkersMat,",
    "              numeric(getNRowStdGenoMultiMarkersMat()))",
    "findIndiceRelatedSample()",
    "kl <- refineKin(cutoff)",
    "saveRDS(list(i = kl$iIndex, j = kl$jIndex, mq = mq, std = std), out)"
  ), script)
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script,
                                 stdout = log, stderr = log))
  testthat::expect_identical(as.integer(st), 0L)
  res <- readRDS(out)
  res$log <- readLines(log, warn = FALSE)
  res
}

# All 1-based (i, j) pairs with i < j, as "i:j" keys.
all_pair_keys <- function(n) {
  p <- utils::combn(n, 2)
  paste(p[1, ], p[2, ], sep = ":")
}

# The pair the loop bound drops, derived by hand from k = C(N,2) - 1:
#   N = 2m   -> i = j = m - 1, folded to (m - 1, m);
#   N = 2m+1 -> (m - 1, 2m), unfolded.
# 1-based, i < j.
dropped_pair_key <- function(n) {
  m <- n %/% 2
  if (n %% 2 == 0) paste(m, m + 1, sep = ":") else paste(m, n, sep = ":")
}

for (n_samples in c(10L, 11L)) {
  label <- paste0("LEDGER #63: the relatedness screen misses one pair at N = ", n_samples)
  test_that(label, {
    skip_unless_level("standard")

    dir <- tempfile()
    dir.create(dir)
    prefix <- file.path(dir, "grmfix")
    write_plink_fixture(prefix, make_geno(40L, n_samples, 63L))
    res <- run_screen(prefix, n_samples)

    expect_identical(res$mq, 40L)             # nothing lost to QC
    n_pairs <- choose(n_samples, 2)

    # The count the code prints is the loop bound itself.
    expect_true(any(grepl(paste0("^totalCombination: ", n_pairs - 1, "$"),
                          res$log)))

    got <- paste(pmin(res$i, res$j), pmax(res$i, res$j), sep = ":")
    expect_identical(anyDuplicated(got), 0L)           # the fold is injective
    expect_true(all(got %in% all_pair_keys(n_samples)))  # and lands in-range

    # ONE pair short, and it is the image of the last k.
    missing <- setdiff(all_pair_keys(n_samples), got)
    expect_identical(length(got), as.integer(n_pairs - 1))
    expect_identical(missing, dropped_pair_key(n_samples))

    # POSITIVE CONTROL that the loop BOUND dropped it, not the relatedness
    # filter: recompute every pair's screening kinship from the same
    # standardised genotypes the worker dots. All C(N,2) of them clear the
    # cutoff, and the missing pair's value sits inside the range of the ones
    # that were returned -- so nothing about the pair itself distinguishes it.
    kin_mat <- crossprod(res$std) / nrow(res$std)
    p_all <- utils::combn(n_samples, 2)
    kin_all <- kin_mat[cbind(p_all[1, ], p_all[2, ])]
    expect_true(all(kin_all > -1e6))
    ij <- as.integer(strsplit(missing, ":")[[1]])
    kin_missing <- kin_mat[ij[1], ij[2]]
    kin_kept <- kin_all[all_pair_keys(n_samples) != missing]
    expect_gte(kin_missing, min(kin_kept))
    expect_lte(kin_missing, max(kin_kept))
  })
}
