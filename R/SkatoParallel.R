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
#
# BATCHED ACROSS REGIONS (step 5, 2026-09-06). Feeding the pool once per
# region left a barrier per region -- the last cell's tail -- which is what
# held the pool at 1.30x against a 1.43x ceiling and made 8 workers no better
# than 4. SAIGE.Region now collects the cells of a whole batch of regions
# (skatoCollect), evaluates them in ONE pool call sorted expensive-first
# (skatoEvaluate) and serves each region back in order (skatoServe); the
# escalated regions' accurate-pass cells go round once more per batch. The
# pieces are separate functions so the driver can interleave the C++ calls
# and the single-variant row copies between them in region order.

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

#' The memory the driver may hold in C++ output lists before it dispatches a
#' batch early (`SAIGE_SKATO_BATCH_MB`, default 512). Step 5 (2026-09-06)
#' evaluates the SKAT-O cells of a whole batch of regions in ONE pool call,
#' which means holding every region's `mainRegionInCPP` output until its
#' answers are back: ~7-10 MB per region and model at N = 50,000 (the N x
#' cells `genoSumMat` dominates), so a 100-region chunk of five models would
#' be several GB. The budget caps the batch below the chunk when needed; a
#' batch of one region is the 2026-09-04 per-region pool.
#' @noRd
skatoBatchBudgetBytes <- function() {
  mb <- suppressWarnings(as.numeric(Sys.getenv("SAIGE_SKATO_BATCH_MB", "512")))
  if (is.na(mb) || mb < 0) mb <- 512
  mb * 1024^2
}

#' COLLECT the SKAT-O cells one `region_postprocess` drive would test.
#'
#' `args` is the argument list of that drive (everything but `skat_fun`). The
#' placeholder answer must be shaped like a real one so the loop's arithmetic
#' (ACAT-V/ACAT-O combination, resultDF construction) runs without error; its
#' VALUES are discarded with the rest of this pass, 0.5 being a legal p-value
#' everywhere one is consumed. Returns the list of (Score, Phi, r.corr,
#' regionTestType) jobs, in call order.
#' @noRd
skatoCollect <- function(args) {
  st <- new.env(parent = emptyenv())
  st$jobs <- list()
  collect <- function(Score, Phi, r.corr, regionTestType) {
    st$jobs[[length(st$jobs) + 1L]] <- list(Score = Score, Phi = Phi, r.corr = r.corr,
                                            regionTestType = regionTestType)
    list(Pvalue_SKATO = 0.5, Pvalue_Burden = 0.5, Pvalue_SKAT = 0.5,
         BETA_Burden = 0, SE_Burden = 1)
  }
  # record_pooled = FALSE: this drive's p-values are the 0.5 placeholders above,
  # so only the SERVE drive may make the LEDGER #127 record.
  invisible(do.call(region_postprocess,
                    c(args, list(skat_fun = collect, record_pooled = FALSE))))
  st$jobs
}

#' Evaluate a batch of collected cells on the pool; answers in job order.
#'
#' One task per cell, `clusterApplyLB` (load-BALANCED: cell cost spans three
#' orders of magnitude, 0.5 ms for a single-marker cell to 540 ms for the
#' worst, and the top 10% of cells hold 48% of the total, so static
#' round-robin would routinely hand one worker most of the work). The cells
#' are dispatched expensive-first -- longest-processing-time-first on the
#' free proxy `ncol(Phi)`, the Davies integration being over Phi's
#' eigenvalues -- so the batch's tail is a small cell, not a large one; the
#' answers are put back in call order before they are served. Ties keep call
#' order (`order` is stable), and the order of EVALUATION cannot reach the
#' bytes anyway: each answer is a pure function of its own cell.
#' @noRd
skatoEvaluate <- function(cl, jobs, dedupe = skatoDedupe(jobs)) {
  if (length(jobs) == 0L) {
    return(list())
  }
  # Identical cells are evaluated once (see skatoDedupe): the pool sees only
  # the distinct jobs and every duplicate is served its twin's answer.
  uniq <- jobs[dedupe$uniq]
  cost <- vapply(uniq, function(j) as.numeric(NCOL(j$Phi)), numeric(1))
  ord <- order(-cost)
  answers <- vector("list", length(uniq))
  answers[ord] <- parallel::clusterApplyLB(cl, uniq[ord], skatoEvalJob)
  skatoDedupeCount(length(jobs), length(uniq))
  answers[dedupe$idx]
}

# ---------------------------------------------------------------------------
# Identical cells evaluated once (2026-09-09).
#
# WHY. The cells of a region are (annotation x max_MAF x weight set), and on
# real data many of them hold exactly the same markers: an annotation union
# `missense;lof;synonymous` equals `missense;lof` in every gene without a
# synonymous marker under the cutoff, a max_MAF band that adds no marker
# repeats the band below it, and so on. In the bundled fixture 4 of the 12
# tested cells per run are such twins (P6: `missense;lof` and
# `missense;lof;synonymous` in both genes). The C++ carry-forward
# (q_maf_for_anno) covers only cutoffs ABOVE an annotation's largest MAF;
# everything else was evaluated again, at ~445 Davies calls a cell.
#
# WHY THIS CANNOT MOVE A BYTE. A cell is served an earlier cell's answer only
# when its (Score, Phi, r.corr, regionTestType) is `identical()` to that
# cell's -- the same bits, attributes included -- and `get_SKAT_pvalue` is a
# pure function of exactly those arguments (the premise the fork pool already
# rests on). Whether the C++ side produces bit-identical inputs for identical
# marker sets is deliberately NOT part of the argument: if it did not, the
# memo would simply never fire, which is a coverage fact and is what the
# printed counts are for. The duplicate ROWS are still emitted, so the
# Cauchy row still combines them exactly as before.
#
# HOW. A cell's key is the hex (`%a`, exact) rendering of its numbers plus
# its shapes; a key hit is confirmed with `identical()` before it counts, so
# the key is only a pre-filter and the contract is R's own equality. Keys are
# compared with `match`, not used as environment names (a 50-marker cell's
# key is ~60 KB, past R's symbol limit).
# ---------------------------------------------------------------------------

#' The exact key of one collected cell.
#' @noRd
skatoJobKey <- function(job) {
  paste(c(job$regionTestType, length(job$Score), dim(job$Score), "|", dim(job$Phi), "|",
          sprintf("%a", c(job$r.corr, job$Score, job$Phi))), collapse = " ")
}

#' Find the identical cells in a batch. Returns `uniq`, the positions of the
#' first occurrence of each distinct cell, and `idx`, for every job the
#' position in `uniq` of the cell it is identical to, so that
#' `answers_for_uniq[idx]` is the answer vector for `jobs`.
#' @noRd
skatoDedupe <- function(jobs) {
  n <- length(jobs)
  if (n == 0L) {
    return(list(uniq = integer(0), idx = integer(0)))
  }
  keys <- vapply(jobs, skatoJobKey, character(1))
  uniq <- integer(0)
  ukeys <- character(0)
  idx <- integer(n)
  for (i in seq_len(n)) {
    hit <- 0L
    for (u in which(ukeys == keys[i])) {
      if (identical(jobs[[i]], jobs[[uniq[u]]])) {
        hit <- u
        break
      }
    }
    if (hit == 0L) {
      uniq <- c(uniq, i)
      ukeys <- c(ukeys, keys[i])
      hit <- length(uniq)
    }
    idx[i] <- hit
  }
  list(uniq = uniq, idx = idx)
}

#' A `skat_fun` that evaluates `kernel` once per distinct cell, for the drives
#' that do not go through the pool (--nThreads=1, and BURDEN): one memo per
#' region_postprocess drive, which is where the twins live.
#' @noRd
skatoMemoKernel <- function(kernel = get_SKAT_pvalue) {
  st <- new.env(parent = emptyenv())
  st$jobs <- list()
  st$keys <- character(0)
  st$answers <- list()
  function(Score, Phi, r.corr, regionTestType) {
    job <- list(Score = Score, Phi = Phi, r.corr = r.corr, regionTestType = regionTestType)
    key <- skatoJobKey(job)
    for (u in which(st$keys == key)) {
      if (identical(job, st$jobs[[u]])) {
        skatoDedupeCount(1L, 0L)
        return(st$answers[[u]])
      }
    }
    ans <- kernel(Score, Phi, r.corr, regionTestType)
    st$jobs[[length(st$jobs) + 1L]] <- job
    st$keys <- c(st$keys, key)
    st$answers[[length(st$answers) + 1L]] <- ans
    skatoDedupeCount(1L, 1L)
    ans
  }
}

# Per-run accounting of the cells asked for and the cells actually evaluated,
# printed once at the end of the run by skatoDedupeSummary(). Always printed,
# so a run where the memo never fired says "0" rather than nothing -- the
# negative control has to be visible too (and a test asserts the positive one).
skato_env$n_cells <- 0L
skato_env$n_evaluated <- 0L

#' @noRd
skatoDedupeReset <- function() {
  skato_env$n_cells <- 0L
  skato_env$n_evaluated <- 0L
  invisible(NULL)
}

#' @noRd
skatoDedupeCount <- function(n_cells, n_evaluated) {
  skato_env$n_cells <- skato_env$n_cells + as.integer(n_cells)
  skato_env$n_evaluated <- skato_env$n_evaluated + as.integer(n_evaluated)
  invisible(NULL)
}

#' @noRd
skatoDedupeSummary <- function() {
  cat("SKAT-O cells: ", skato_env$n_cells, " asked for, ", skato_env$n_evaluated,
      " evaluated, ", skato_env$n_cells - skato_env$n_evaluated,
      " served from an identical cell.\n", sep = "")
  invisible(c(cells = skato_env$n_cells, evaluated = skato_env$n_evaluated))
}

#' SERVE a batch's answers back into the same `region_postprocess` drive.
#'
#' The cell loop's control flow does not depend on any p-value, so this drive
#' asks for exactly the cells the collecting drive found, in the same order.
#' The counter checks are not decoration: if that ever stopped being true it
#' would otherwise silently pair each cell with another cell's p-value.
#' @noRd
skatoServe <- function(args, answers) {
  st <- new.env(parent = emptyenv())
  st$i <- 0L
  serve <- function(Score, Phi, r.corr, regionTestType) {
    st$i <- st$i + 1L
    if (st$i > length(answers)) {
      stop("skatoServe: the serving pass asked for more SKAT-O cells (", st$i,
           ") than the collecting pass found (", length(answers),
           "). The cell loop's control flow must not depend on a p-value.\n")
    }
    answers[[st$i]]
  }
  out <- do.call(region_postprocess, c(args, list(skat_fun = serve)))
  if (st$i != length(answers)) {
    stop("skatoServe: the serving pass asked for ", st$i,
         " SKAT-O cells but the collecting pass found ", length(answers), ".\n")
  }
  out
}

#' Drive `region_postprocess` with its SKAT-O cells evaluated in parallel,
#' for ONE region: collect, evaluate, serve. The driver batches regions
#' itself (step 5); this composition is kept for a single region and for the
#' unit tests of the three pieces. `...` is passed through verbatim, so this
#' stays correct as `region_postprocess`'s signature grows.
#' @noRd
region_postprocess_parallel <- function(..., nThreads) {
  if (!canParallelSkato(nThreads)) {
    return(region_postprocess(...))
  }
  args <- list(...)
  jobs <- skatoCollect(args)
  # Nothing to test (a region with no polymorphic cell), or so little that the
  # round trip would cost more than it saves: fall through to the serial path,
  # which also re-runs the loop only once.
  if (length(jobs) < 2L) {
    return(region_postprocess(...))
  }
  skatoServe(args, skatoEvaluate(startSkatoCluster(nThreads), jobs))
}
