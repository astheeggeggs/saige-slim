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
  fake_post <- function(tag, skat_fun = SAIGE:::get_SKAT_pvalue) {
    a <- skat_fun(c(1, 2), diag(2), 0, "SKAT-O")
    b <- skat_fun(3, matrix(4), 0, "SKAT-O")
    list(pval.Region = list(tag = tag, a = a, b = b))
  }
  testthat::local_mocked_bindings(region_postprocess = fake_post, .package = "SAIGE")

  jobs <- SAIGE:::skatoCollect(list("t1"))
  expect_length(jobs, 2L)
  expect_identical(jobs[[1]]$Score, c(1, 2))
  expect_identical(jobs[[2]]$Phi, matrix(4))

  answers <- list(list(p = "first"), list(p = "second"))
  out <- SAIGE:::skatoServe(list("t1"), answers)
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
