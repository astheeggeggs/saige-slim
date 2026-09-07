# Characterization of LEDGER #103: on the SKAT-O path BETA_Burden and SE_Burden
# do not estimate the burden effect. They estimate kappa_w times it, where
# kappa_w = 1'Phi 1 / tr(Phi) is the LD inflation of the weighted score
# covariance.
#
# The score for the burden coefficient b is S_w = 1'Score and its information is
# B_w = 1'Phi 1, so the estimator is S_w / B_w with SE = 1 / sqrt(B_w). SAIGE
# divides by A_w = tr(Phi) instead, giving E[BETA_Burden] = b * kappa_w.
#
# This is NOT an internal inconsistency: Var(S_w / A_w) = B_w / A_w^2, so the
# reported SE_Burden = sqrt(B_w) / A_w is the CORRECT standard error of the
# statistic actually reported, and BETA/SE reproduces the printed p exactly.
# The defect is the estimand, which is why nothing downstream can detect it.
# (The planted-effect and cross-mode reconciliations are in the ledger row.)
#
# These tests assert what the code does TODAY. When #103 is fixed, invert the
# assertions marked INVERT below: BETA_Burden becomes sum(Score) / sum(Phi) and
# SE_Burden becomes 1 / sqrt(sum(Phi)). The others hold either way.

# A small PSD Phi with a large, unambiguous off-diagonal mass: AR(1) correlation
# 0.5 on five markers, so sum(Phi) is ~2x its trace and the two candidate
# denominators cannot be confused for one another at any plausible tolerance.
burden_trace_case <- function() {
  m <- 5
  corMat <- 0.5^abs(outer(seq_len(m), seq_len(m), "-"))
  sdVec <- c(1.3, 0.8, 2.1, 1.1, 0.6)
  list(
    Score = c(9, -4, 12, 3, -2),
    Phi = (sdVec %o% sdVec) * corMat * 100,
    r.corr = SKAT:::SKAT_Check_Method(method = "optimal.adj", r.corr = 0)$r.corr
  )
}

test_that("r.corr = 0 expands to the 7-point SKAT-O grid these tests rely on", {
  expect_equal(burden_trace_case()$r.corr, c(0, 0.01, 0.04, 0.09, 0.25, 0.5, 1))
})

test_that("BETA_Burden divides by the trace of Phi, not by its full sum", {
  tc <- burden_trace_case()
  out <- SAIGE:::get_SKAT_pvalue(tc$Score, tc$Phi, tc$r.corr, "SKAT-O")

  traceOnly <- sum(diag(tc$Phi))
  fullSum <- sum(tc$Phi)
  # the fixture is only meaningful while the two denominators are far apart
  expect_gt(fullSum / traceOnly, 2)

  # INVERT ON FIX. What it does today, to the last bit:
  expect_identical(out$BETA_Burden, sum(tc$Score) / traceOnly)
  # INVERT ON FIX. And it is NOT the estimator of the burden coefficient:
  expect_false(isTRUE(all.equal(out$BETA_Burden, sum(tc$Score) / fullSum)))
  # INVERT ON FIX. The bias is exactly kappa_w, which is what #103 asserts;
  # after the fix this ratio is 1, not kappa_w.
  expect_equal(out$BETA_Burden / (sum(tc$Score) / fullSum), fullSum / traceOnly,
               tolerance = 1e-12)
})

test_that("the Pvalue_Burden printed beside it uses the FULL sum of Phi", {
  # The same row's p-value comes from SKAT's rho = 1 grid element, which uses
  # sum(Phi) -- the correct information for the burden coefficient. Inverting it
  # recovers a variance matching the full sum and not the trace. This is what
  # makes SE_Burden right for the beta that is reported: the p-value supplies
  # B_w, the beta supplies A_w, and the back-derivation combines them into
  # sqrt(B_w)/A_w. (The residual is Davies at rho = 0.999 versus a normal
  # approximation, not the choice of denominator.)
  tc <- burden_trace_case()
  out <- SAIGE:::get_SKAT_pvalue(tc$Score, tc$Phi, tc$r.corr, "SKAT-O")

  impliedVar <- sum(tc$Score)^2 / qchisq(out$Pvalue_Burden, df = 1, lower.tail = FALSE)
  expect_equal(impliedVar / sum(tc$Phi), 1, tolerance = 0.01)
  expect_gt(impliedVar / sum(diag(tc$Phi)), 2)
})

test_that("SE_Burden is the correct standard error OF THE REPORTED beta", {
  # sqrt(B_w)/A_w, not 1/sqrt(B_w): SAIGE's SE is right for SAIGE's beta, which
  # is precisely why #103 is an estimand defect and not an arithmetic one.
  # INVERT ON FIX: SE_Burden becomes 1/sqrt(sum(Phi)).
  # Toleranced at 1e-3: SE is back-derived through qnorm from SKAT's Davies p at
  # rho = 0.999, which is not exactly the normal-approximation p.
  tc <- burden_trace_case()
  out <- SAIGE:::get_SKAT_pvalue(tc$Score, tc$Phi, tc$r.corr, "SKAT-O")

  expect_equal(out$SE_Burden, sqrt(sum(tc$Phi)) / sum(diag(tc$Phi)), tolerance = 1e-3)
  expect_false(isTRUE(all.equal(out$SE_Burden, 1 / sqrt(sum(tc$Phi)), tolerance = 1e-3)))
})

test_that("the estimand defect is invisible to any Z-based check on the row", {
  # SE_Burden is back-derived as abs(BETA_Burden / qnorm(p / 2)), so BETA and SE
  # carry the same scale error and BETA/SE reproduces the printed p exactly.
  # This is why no golden invariant or meta-analysis Z can detect #103, and why
  # the fix moves BETA_Burden and SE_Burden while every p-value stays identical.
  tc <- burden_trace_case()
  out <- SAIGE:::get_SKAT_pvalue(tc$Score, tc$Phi, tc$r.corr, "SKAT-O")

  # Toleranced, not bit-identical: the two sides reach the same quantity by
  # different arithmetic, so the last ULP is the platform's, not the estimator's.
  # This assertion holds both before and after the #103 fix -- do not invert it.
  expect_equal(
    abs(out$BETA_Burden / out$SE_Burden),
    abs(qnorm(out$Pvalue_Burden / 2)),
    tolerance = 1e-12
  )
})

test_that("a single-marker cell cannot move: its trace IS its full sum", {
  # Why the #103 fix is predicted to leave every single-effective-marker row
  # byte-identical, and only multi-marker cells re-baseline.
  tc <- burden_trace_case()
  Phi1 <- tc$Phi[1, 1, drop = FALSE]
  expect_identical(sum(diag(Phi1)), sum(Phi1))

  out <- SAIGE:::get_SKAT_pvalue(tc$Score[1], Phi1, tc$r.corr, "SKAT-O")
  expect_identical(out$BETA_Burden, tc$Score[1] / sum(Phi1))
})
