# Parallel evaluation of a region's SKAT-O p-values.
#
# WHY THIS EXISTS. `Rprof` at N = 50,000 over 210 genes puts 53.8% of a group
# run in `region_postprocess` and 44.9% in `SKAT_davies` under `integrate()` --
# R adaptive quadrature over the SKAT package's compiled Davies routine, ~445
# calls per (annotation x max_MAF x weight set) cell at ~94 us each. The C++
# region loop is the other 43.9%, so threading THAT caps at 1.49x on 4 threads
# for SKAT-O. See `harness/bench/NOTES.md`, "Where the OTHER 90% of a gene test
# goes".
#
# WHY NOT THE OTHER LEVERS, both measured and both rejected (same doc):
#   * rewriting in C++ buys ~1.8% of the run -- `SKAT_davies` is already a thin
#     wrapper around `.C("qfc")`, so the hot code is compiled already;
#   * loosening the accuracy SKAT asks of Davies (`acc = 1e-6`) is worth 24-30%
#     and moves genes with p < 1e-4 by up to 80% relative.
#
# WHY A PERSISTENT CLUSTER RATHER THAN `mclapply`. Measured on the real
# workload: forking per region costs 55-90 ms against ~350 ms of work and is a
# net SLOWDOWN at 8 workers (0.91x). A cluster forked ONCE costs 8-20 ms a
# region and gives ~3x on the parallel part.
#
# WHY THE RESULT IS BIT-IDENTICAL, not merely close. Each cell's SKAT-O p-value
# is a pure function of that cell's (Score, Phi, r.corr) -- no cross-cell
# reduction, no shared accumulator, no RNG. The workers run the SAME
# `get_SKAT_pvalue` on the SAME inputs, so this is one code path re-run over a
# partition, which is the narrow case where byte-identity is the right gate
# (cf. CLAUDE.md's standing rule against comparing two different code paths).
#
# HOW. `region_postprocess`'s cell loop is a 175-line triple nest with the two
# SKAT calls buried in conditional branches, and its control flow does NOT
# depend on any p-value it computes. So rather than restructure it, it is driven
# TWICE through its injected `skat_fun`: once to collect every (Score, Phi) it
# would test, then the collected set is evaluated in parallel, then again to
# serve the answers back in call order. The second drive costs the loop's
# non-SKAT work a second time; that is the price of not rewriting the most
# intricate function in the package, and it is measured in NOTES.md.

# One worker pool for the whole run. Created on first use, reused by every
# region, and torn down by stopSkatoCluster().
skato_env <- new.env(parent = emptyenv())
skato_env$cl <- NULL
skato_env$n <- 0L

#' Start (or reuse) the SKAT-O worker pool.
#' @noRd
startSkatoCluster <- function(nThreads) {
  if (!is.null(skato_env$cl) && identical(skato_env$n, as.integer(nThreads))) {
    return(invisible(skato_env$cl))
  }
  stopSkatoCluster()
  # A FORK cluster: the workers inherit SKAT and the loaded model, so nothing
  # has to be exported, and the tasks carry only (Score, Phi) -- a few KB for a
  # typical cell. PSOCK measured marginally faster here but needs the package
  # loaded and every object shipped, so fork is the simpler contract. It is
  # also why this is Unix-only; see the caller's guard.
  skato_env$cl <- parallel::makeForkCluster(as.integer(nThreads))
  skato_env$n <- as.integer(nThreads)
  # Printed exactly once, and only once the pool is actually built -- i.e. only
  # when a region really did take the parallel path. A test that asserts
  # byte-identity without asserting THIS would pass just as happily if the
  # option never arrived and every region fell through to the serial branch,
  # which is the LEDGER #30 shape.
  cat("SKAT-O post-processing will use ", skato_env$n, " workers\n", sep = "")
  invisible(skato_env$cl)
}

#' Tear the pool down. Safe to call when there is none.
#' @noRd
stopSkatoCluster <- function() {
  if (!is.null(skato_env$cl)) {
    try(parallel::stopCluster(skato_env$cl), silent = TRUE)
  }
  skato_env$cl <- NULL
  skato_env$n <- 0L
  invisible(NULL)
}

#' Evaluate one collected SKAT-O cell.
#'
#' Defined at TOP LEVEL on purpose, and this is not a style point. R serializes
#' a closure together with its enclosing environment, so a worker function
#' defined inside `region_postprocess_parallel` would drag that frame's `...`
#' -- the whole `outList`, `VarMat` and all -- to a worker on EVERY task
#' dispatch. Measured: that cost ~13 ms per task against a ~40 ms task and made
#' 4 workers SLOWER than serial at N = 50,000. A top-level function's
#' environment is the package namespace, which the workers already have, so
#' nothing but `job` is shipped.
#' @noRd
skatoEvalJob <- function(job) {
  get_SKAT_pvalue(job$Score, job$Phi, job$r.corr, job$regionTestType)
}

#' TRUE when a parallel SKAT-O pass is possible and asked for.
#' @noRd
canParallelSkato <- function(nThreads) {
  is.numeric(nThreads) && length(nThreads) == 1L && !is.na(nThreads) &&
    nThreads > 1 && .Platform$OS.type == "unix"
}

#' Drive `region_postprocess` with its SKAT-O cells evaluated in parallel.
#'
#' `...` is passed through verbatim, so this stays correct as
#' `region_postprocess`'s signature grows.
#' @noRd
region_postprocess_parallel <- function(..., nThreads) {
  if (!canParallelSkato(nThreads)) {
    return(region_postprocess(...))
  }

  # --- pass 1: collect -----------------------------------------------------
  # The placeholder must be shaped like a real answer so the loop's arithmetic
  # (ACAT-V/ACAT-O combination, resultDF construction) runs without error. Its
  # VALUES are discarded with the rest of this pass; 0.5 is used because it is
  # a legal p-value everywhere one is consumed.
  st <- new.env(parent = emptyenv())
  st$jobs <- list()
  collect <- function(Score, Phi, r.corr, regionTestType) {
    st$jobs[[length(st$jobs) + 1L]] <- list(Score = Score, Phi = Phi, r.corr = r.corr,
                                            regionTestType = regionTestType)
    list(Pvalue_SKATO = 0.5, Pvalue_Burden = 0.5, Pvalue_SKAT = 0.5,
         BETA_Burden = 0, SE_Burden = 1)
  }
  invisible(region_postprocess(..., skat_fun = collect))

  # Nothing to test (a region with no polymorphic cell), or so little that the
  # round trip would cost more than it saves: fall through to the serial path,
  # which also re-runs the loop only once.
  if (length(st$jobs) < 2L) {
    return(region_postprocess(...))
  }

  # --- pass 2: evaluate in parallel ---------------------------------------
  cl <- startSkatoCluster(nThreads)
  # Load-BALANCED: cell cost spans three orders of magnitude (0.5 ms for a
  # single-marker cell to 540 ms for the worst), and the top 10% of cells hold
  # 48% of the total, so static round-robin scheduling would routinely hand one
  # worker most of the work.
  answers <- parallel::clusterApplyLB(cl, st$jobs, skatoEvalJob)

  # --- pass 3: serve -------------------------------------------------------
  # The cell loop's control flow does not depend on any p-value, so this drive
  # asks for exactly the same cells in exactly the same order. The counter
  # check is not decoration: if that ever stopped being true it would otherwise
  # silently pair each cell with another cell's p-value.
  st$i <- 0L
  serve <- function(Score, Phi, r.corr, regionTestType) {
    st$i <- st$i + 1L
    if (st$i > length(answers)) {
      stop("region_postprocess_parallel: the serving pass asked for more SKAT-O ",
           "cells (", st$i, ") than the collecting pass found (", length(answers),
           "). The cell loop's control flow must not depend on a p-value.\n")
    }
    answers[[st$i]]
  }
  out <- region_postprocess(..., skat_fun = serve)
  if (st$i != length(answers)) {
    stop("region_postprocess_parallel: the serving pass asked for ", st$i,
         " SKAT-O cells but the collecting pass found ", length(answers), ".\n")
  }
  out
}
