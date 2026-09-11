# The three pieces of the batched SKAT-O pool (step 5, 2026-09-06;
# R/SkatoParallel.R): skatoCollect gathers the cells one region_postprocess
# drive would test, skatoEvaluate answers a batch on the pool expensive-first
# and puts the answers back in job order, skatoServe hands them back into the
# same drive and refuses a count mismatch. Byte-identity of a whole run is
# test-step2-skato-parallel.R; this pins the pieces on small cells built by
# hand, against the serial kernel.
#
# Tier: fast (a 2-worker fork pool for a few milliseconds; unix only).

skip_on_os("windows")

# A well-posed SKAT-O cell of m markers: Phi symmetric positive definite.
ub_cell <- function(m, seed) {
  set.seed(seed)
  A <- matrix(rnorm(m * m), m, m)
  Phi <- crossprod(A) + diag(m) * 0.5
  list(Score = rnorm(m, sd = 0.5), Phi = Phi, r.corr = 0, regionTestType = "SKAT-O")
}

test_that("skatoEvaluate returns the serial answers in job order whatever the sort", {
  # sizes chosen so the expensive-first order is NOT the call order
  sizes <- c(2L, 7L, 1L, 5L, 7L, 3L)
  jobs <- Map(ub_cell, sizes, seq_along(sizes))
  serial <- lapply(jobs, SAIGE:::skatoEvalJob)
  cl <- parallel::makeForkCluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  batched <- SAIGE:::skatoEvaluate(cl, jobs)
  expect_identical(batched, serial)
  # an empty batch is an empty answer list, not a pool call
  expect_identical(SAIGE:::skatoEvaluate(cl, list()), list())
})

test_that("skatoCollect gathers the drive's cells and skatoServe checks the count", {
  # a stand-in for region_postprocess that asks for two cells and returns
  # what it was told, so the served answers are visible in the result
  # `record_pooled` is captured rather than ignored: LEDGER #127 records the
  # pooled row's provenance from inside region_postprocess, and the parallel
  # path drives that body TWICE, so only the SERVE drive may record. A
  # collecting drive that recorded would double-count every region silently.
  seen <- new.env(parent = emptyenv())
  seen$record <- list()
  fake_post <- function(tag, skat_fun = SAIGE:::get_SKAT_pvalue, record_pooled = TRUE) {
    seen$record[[length(seen$record) + 1L]] <- record_pooled
    a <- skat_fun(c(1, 2), diag(2), 0, "SKAT-O")
    b <- skat_fun(3, matrix(4), 0, "SKAT-O")
    list(pval.Region = list(tag = tag, a = a, b = b))
  }
  testthat::local_mocked_bindings(region_postprocess = fake_post, .package = "SAIGE")

  jobs <- SAIGE:::skatoCollect(list("t1"))
  expect_identical(seen$record, list(FALSE))
  expect_length(jobs, 2L)
  expect_identical(jobs[[1]]$Score, c(1, 2))
  expect_identical(jobs[[2]]$Phi, matrix(4))

  answers <- list(list(p = "first"), list(p = "second"))
  out <- SAIGE:::skatoServe(list("t1"), answers)
  expect_identical(seen$record, list(FALSE, TRUE))   # collect, then serve
  expect_identical(out$pval.Region$tag, "t1")
  expect_identical(out$pval.Region$a, answers[[1]])
  expect_identical(out$pval.Region$b, answers[[2]])

  expect_error(SAIGE:::skatoServe(list("t1"), answers[1]), "asked for more SKAT-O cells")
  expect_error(SAIGE:::skatoServe(list("t1"), c(answers, answers[1])),
               "but the collecting pass found 3")
})

test_that("the batch memory budget reads the environment and rejects nonsense", {
  withr::with_envvar(c(SAIGE_SKATO_BATCH_MB = "0"), {
    expect_identical(SAIGE:::skatoBatchBudgetBytes(), 0)
  })
  withr::with_envvar(c(SAIGE_SKATO_BATCH_MB = "2"), {
    expect_identical(SAIGE:::skatoBatchBudgetBytes(), 2 * 1024^2)
  })
  withr::with_envvar(c(SAIGE_SKATO_BATCH_MB = "banana"), {
    expect_identical(SAIGE:::skatoBatchBudgetBytes(), 512 * 1024^2)
  })
})

# ---------------------------------------------------------------------------
# Identical cells evaluated once (2026-09-09): skatoDedupe / skatoMemoKernel /
# the counters. A cell is served an earlier cell's answer only when its inputs
# are identical(); both controls are pinned -- a one-ULP perturbation and an
# attribute-only difference must each stay distinct.
# ---------------------------------------------------------------------------

test_that("skatoDedupe pairs identical cells and refuses near-identical ones", {
  A <- ub_cell(3L, 11)
  B <- ub_cell(3L, 12)
  C <- ub_cell(2L, 13)
  A_ulp <- A
  # (1 + 2 eps) moves any double by at least one ULP; (1 + eps) can round back
  A_ulp$Score[1] <- A_ulp$Score[1] * (1 + 2 * .Machine$double.eps)
  A_attr <- A
  dimnames(A_attr$Phi) <- list(letters[1:3], letters[1:3])         # same numbers, extra attribute
  jobs <- list(A, B, A, C, B, A_ulp, A_attr, A)
  dd <- SAIGE:::skatoDedupe(jobs)
  expect_identical(dd$uniq, c(1L, 2L, 4L, 6L, 7L))
  expect_identical(dd$idx, c(1L, 2L, 1L, 3L, 2L, 4L, 5L, 1L))
  # the fan-out reproduces the job list exactly
  expect_identical(jobs[dd$uniq][dd$idx], jobs)
  expect_identical(SAIGE:::skatoDedupe(list()), list(uniq = integer(0), idx = integer(0)))
  # the key itself is exact: the ULP twin and the attribute twin share A's numbers
  expect_identical(SAIGE:::skatoJobKey(A_attr), SAIGE:::skatoJobKey(A))
  expect_false(identical(SAIGE:::skatoJobKey(A_ulp), SAIGE:::skatoJobKey(A)))
})

test_that("skatoEvaluate serves duplicates their twin's answer and counts them", {
  jobs <- list(ub_cell(4L, 21), ub_cell(2L, 22), ub_cell(4L, 21), ub_cell(6L, 23), ub_cell(2L, 22))
  serial <- lapply(jobs, SAIGE:::skatoEvalJob)
  cl <- parallel::makeForkCluster(2L)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  SAIGE:::skatoDedupeReset()
  batched <- SAIGE:::skatoEvaluate(cl, jobs)
  expect_identical(batched, serial)
  counts <- SAIGE:::skatoDedupeSummary()
  expect_identical(unname(counts), c(5L, 3L))
  expect_output(SAIGE:::skatoDedupeSummary(),
                "SKAT-O cells: 5 asked for, 3 evaluated, 2 served from an identical cell.",
                fixed = TRUE)
})

test_that("skatoMemoKernel calls the kernel once per distinct cell", {
  ct <- new.env()
  ct$calls <- 0L
  counting <- function(Score, Phi, r.corr, regionTestType) {
    ct$calls <- ct$calls + 1L
    SAIGE:::get_SKAT_pvalue(Score, Phi, r.corr, regionTestType)
  }
  A <- ub_cell(3L, 31)
  B <- ub_cell(5L, 32)
  A_ulp <- A
  A_ulp$Phi[2, 2] <- A_ulp$Phi[2, 2] * (1 + 2 * .Machine$double.eps)
  ask <- function(f, j) f(j$Score, j$Phi, j$r.corr, j$regionTestType)
  SAIGE:::skatoDedupeReset()
  memo <- SAIGE:::skatoMemoKernel(counting)
  got <- lapply(list(A, B, A, A_ulp, B, A), function(j) ask(memo, j))
  expect_identical(ct$calls, 3L)              # A, B, A_ulp, and no more
  want <- lapply(list(A, B, A, A_ulp, B, A), function(j) ask(SAIGE:::get_SKAT_pvalue, j))
  expect_identical(got, want)
  expect_identical(unname(SAIGE:::skatoDedupeSummary()), c(6L, 3L))
  # a fresh memo remembers nothing (one per region_postprocess drive)
  memo2 <- SAIGE:::skatoMemoKernel(counting)
  ask(memo2, A)
  expect_identical(ct$calls, 4L)
})
