# Base-R builtins whose SEMANTICS -- not merely their result -- change with the
# type, length or shape of the argument. These are the idioms behind LEDGER
# #33-#37; the tests pin the language behaviour itself, so they stay meaningful
# after those call sites are fixed and they cost nothing to run.
#
# See the "Type-polymorphic idiom inventory" section of LEDGER.md.

test_that("diag() dispatches on the type and shape of its argument", {
  # the intended use: extract the diagonal of a square matrix
  m <- matrix(c(1, 2, 3, 4), 2, 2)
  expect_identical(diag(m), c(1, 4))

  # a q x 1 matrix silently yields only the FIRST element, not the column.
  # This is why R/SAIGE_SPATest_Region_Func.R:250 wrapping a vector in
  # as.matrix() is a defence only when the argument is already a matrix.
  expect_identical(diag(as.matrix(c(1, 2, 3))), 1)

  # a length > 1 vector CONSTRUCTS a diagonal matrix instead
  expect_identical(dim(diag(c(0.2, 0.3))), c(2L, 2L))
  expect_identical(sum(diag(c(0.2, 0.3))), 0.5)

  # a length-1 numeric is read as a DIMENSION: diag(x) is a trunc(x)-square
  # identity matrix. Below 1 that is 0 x 0, so sum() is 0 -- the mechanism
  # behind LEDGER #34's Inf beta.
  expect_identical(dim(diag(0.37)), c(0L, 0L))
  expect_identical(sum(diag(0.37)), 0)
  expect_identical(dim(diag(3.7)), c(3L, 3L))
  expect_identical(sum(diag(3.7)), 3)

  # a 1 x 1 matrix -- what drop = FALSE preserves -- behaves correctly
  expect_identical(diag(matrix(0.37, 1, 1)), 0.37)
  expect_identical(sum(diag(matrix(0.37, 1, 1))), 0.37)
})

test_that("dropping a single-element matrix subset is what feeds diag() a scalar", {
  m <- matrix(c(0.42, 0.1, 0.1, 0.7), 2, 2)
  pos <- 1L
  expect_null(dim(m[pos, pos]))                       # LEDGER #34, line 455
  expect_identical(dim(m[pos, pos, drop = FALSE]), c(1L, 1L))
  # the two disagree only once diag() is applied
  expect_identical(sum(diag(m[pos, pos])), 0)
  expect_identical(sum(diag(m[pos, pos, drop = FALSE])), 0.42)
})

test_that("matrix - vector recycles column-major, not down the diagonal", {
  # LEDGER #33: subtracting the extracted diagonal from a matrix applies the
  # i-th diagonal entry across the whole of row i, and the result is not
  # symmetric. Only the diagonal itself survives intact, which is why
  # BETA_Burden_cond is unaffected but the SKAT covariance is not.
  P <- matrix(0, 3, 3)
  V <- matrix(1:9, 3, 3)
  got <- P - diag(V)
  expect_identical(got, matrix(rep(c(-1, -5, -9), 3), 3, 3))
  expect_false(identical(got, P - V))
  expect_identical(diag(got), diag(P - V))            # diagonals do agree
  # ... and the two coincide exactly at q == 1
  P1 <- matrix(2, 1, 1)
  V1 <- matrix(0.5, 1, 1)
  expect_identical(as.numeric(P1 - diag(V1)), as.numeric(P1 - V1))
})

test_that("sample() reads a length-1 first argument as a range, not as the set", {
  # LEDGER inventory: which() returning exactly one index turns
  # sample(idx, 1) into a uniform draw from 1:idx.
  set.seed(1)
  draws <- vapply(1:200, function(i) sample(5000L, size = 1), numeric(1))
  expect_false(all(draws == 5000L))
  expect_true(any(draws < 5000L))
  # the safe idiom: index into the vector rather than passing it as the set
  x <- 5000L
  expect_identical(x[sample.int(length(x))], x)
  # length > 1 behaves as expected, which is what hides the bug in testing
  y <- c(11L, 22L)
  expect_setequal(sample(y, size = 2), y)
})

test_that("1:n runs backwards at n == 0 where seq_len() is empty", {
  # LEDGER #35, the categorical variance-ratio loop with one MAC category,
  # and LEDGER #37.
  expect_identical(1:0, c(1L, 0L))
  expect_length(1:0, 2L)
  expect_length(seq_len(0), 0L)
  # seq(1:n) -- R/SAIGE_fitGLMM_fast.R:1964 -- is doubly wrong at n == 0
  expect_identical(seq(1:0), 1:2)
  # and seq(1, 0, by = 1) errors before any hand-written guard can fire
  expect_error(seq(1, 0, by = 1), "wrong sign in 'by' argument")
})

test_that("rbind recycles a short row instead of rejecting it", {
  # LEDGER #38: a --X_PARregion token with no "-" yields length 1.
  m <- rbind(matrix(c(10, 20), 1, 2), as.numeric(unlist(strsplit("5", "-"))))
  expect_identical(dim(m), c(2L, 2L))
  expect_identical(as.numeric(m[2, ]), c(5, 5))
})
