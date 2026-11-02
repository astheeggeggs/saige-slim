# Step-0 fixtures and the in-process screen/refine runner, shared by
# test-step0-blocked-screen.R and test-step0-refine-pairmajor.R (moved out of
# the former on 2026-11-02; testthat runs files in separate processes, so a
# test file cannot source another).

write_plink_fixture_b <- function(prefix, geno) {
  m <- nrow(geno)
  n <- ncol(geno)
  writeLines(sprintf("1\tsnp%d\t0\t%d\tA\tT", seq_len(m), seq_len(m)),
             paste0(prefix, ".bim"))
  writeLines(sprintf("f%d\ts%d\t0\t0\t0\t-9", seq_len(n), seq_len(n)),
             paste0(prefix, ".fam"))
  code <- c(3L, 2L, 0L)
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

# Haplotype fixture: founders, then `npo` parent-offspring pairs and `ndup`
# duplicates appended. Every marker polymorphic.
make_related_geno <- function(m, nfounder, npo, ndup, seed) {
  set.seed(seed)
  p <- stats::runif(m, 0.05, 0.5)
  h1 <- matrix(stats::rbinom(m * nfounder, 1, p), m, nfounder)
  h2 <- matrix(stats::rbinom(m * nfounder, 1, p), m, nfounder)
  par <- seq_len(npo)
  pick <- matrix(stats::runif(m * npo) < 0.5, m, npo)
  kid1 <- ifelse(pick, h1[, par, drop = FALSE], h2[, par, drop = FALSE])
  kid2 <- matrix(stats::rbinom(m * npo, 1, p), m, npo)
  dup <- seq_len(ndup) + npo
  g <- cbind(h1 + h2, kid1 + kid2, (h1 + h2)[, dup, drop = FALSE])
  for (i in seq_len(m)) {
    if (sum(g[i, ]) == 0) g[i, 1] <- 1L
    if (sum(g[i, ]) == 2 * ncol(g)) g[i, 1] <- 1L
  }
  list(g = g,
       planted = c(paste(par, nfounder + par, sep = ":"),
                   paste(dup, nfounder + npo + seq_len(ndup), sep = ":")))
}

# One child process per screen: setgeno and both screens write C++ globals.
# Returns the screened pairs as sorted "i:j" keys (1-based, as refineKin
# returns them) plus the log.
run_screen_b <- function(prefix, n, method, block = 2048L, cutoff = -1e6,
                         threads = 1L, refine = c("pairmajor", "markermajor"),
                         refine_block = NULL) {
  refine <- match.arg(refine)
  out <- tempfile(fileext = ".rds")
  log <- tempfile(fileext = ".log")
  script <- tempfile(fileext = ".R")
  call <- if (method == "blocked") {
    sprintf("findIndiceRelatedSampleBlocked(%dL)", as.integer(block))
  } else {
    "findIndiceRelatedSample()"
  }
  writeLines(c(
    "suppressPackageStartupMessages(library(SAIGE))",
    sprintf('prefix <- "%s"; out <- "%s"; n <- %d; cutoff <- %g',
            prefix, out, n, cutoff),
    sprintf("RcppParallel::setThreadOptions(numThreads = %d)", as.integer(threads)),
    "setminMAFforGRM(0.01); setmaxMissingRateforGRM(0.15)",
    'setgeno(paste0(prefix, ".bed"), paste0(prefix, ".bim"),',
    '        paste0(prefix, ".fam"), seq_len(n), rep(TRUE, n), 2, FALSE)',
    "mq <- length(getAlleleFreqVec())",
    "setRelatednessCutoff(cutoff)",
    "setSubMarkerIndex(seq_len(mq) - 1L)",
    "Get_MultiMarkersBySample_StdGeno_Mat()",
    call,
    # refineKin at -1e30 returns every screened pair; refineKin itself is
    # unchanged and shared by both screens.
    # the refine switch and block length (LEDGER #185) are read from the
    # environment at the refineKin call
    sprintf('Sys.setenv(SAIGE_REFINE_PAIRMAJOR = "%s")', if (refine == "pairmajor") "1" else "0"),
    if (is.null(refine_block)) "Sys.unsetenv(\"SAIGE_REFINE_BLOCK\")" else
      sprintf('Sys.setenv(SAIGE_REFINE_BLOCK = "%d")', as.integer(refine_block)),
    "kl <- refineKin(-1e30)",
    "saveRDS(list(i = kl$iIndex, j = kl$jIndex, k = kl$kinValue), out)"
  ), script)
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script,
                                 stdout = log, stderr = log))
  testthat::expect_identical(as.integer(st), 0L)
  res <- readRDS(out)
  o <- order(res$i, res$j)
  list(i = res$i[o], j = res$j[o], k = res$k[o],
       keys = paste(res$i[o], res$j[o], sep = ":"),
       log = readLines(log, warn = FALSE))
}
