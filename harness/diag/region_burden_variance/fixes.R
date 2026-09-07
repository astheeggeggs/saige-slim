# The candidate fixes for LEDGER #103/#104/#105 as functions over the installed
# SAIGE namespace, so the pipeline driver (drv.R, via SAIGE_FIX) and the
# offline analysis (sim_analyze.R) run one and the same code. Nothing here is
# installed unless a caller does so.
diag_fixes <- function(ns = asNamespace("SAIGE")) {
  orig_skat <- get("get_SKAT_pvalue", envir = ns)
  orig_phi <- get("get_newPhi_scaleFactor_traitType", envir = ns)
  orig_adj <- get("SPA_ER_kernel_related_Phiadj_fast_new", envir = ns)
  spa_pval <- get("SPA_pval", envir = ns)
  # Replace the one top-level statement of f that assigns `lhs`, checking that
  # its right-hand side is the expected one (deparse formatting is not relied on).
  edit_stmt <- function(f, lhs, expect_rhs, new_stmt) {
    b <- body(f)
    hit <- vapply(seq_along(b), function(i) {
      s <- b[[i]]
      is.call(s) && identical(s[[1]], as.name("<-")) && identical(s[[2]], as.name(lhs)) &&
        grepl(expect_rhs, paste(deparse(s[[3]], width.cutoff = 500L), collapse = " "), fixed = TRUE)
    }, logical(1))
    stopifnot(sum(hit) == 1)
    b[[which(hit)]] <- new_stmt
    body(f) <- b
    environment(f) <- ns
    f
  }
  # LEDGER #103: the burden information is the full sum of Phi, not its trace.
  skat103 <- edit_stmt(orig_skat, "varSum_Phi", "sum(diag(Phi))",
                       quote(varSum_Phi <- sum(Phi)))
  # Fix A (#104/#105): standardise the burden statistic by the mixed-model
  # burden variance 1'Phi 1 before the |z| < 2 test and the saddlepoint -- the
  # device the C++ marker SPA uses (q = Tstat/sqrt(var1/var2) + m1). Below the
  # cutoff VarQ_2 then equals 1'Phi 1 and the r rescale is a no-op.
  phi_a <- function(q.sum, mu.a, g.sum, p.new, Score, Phi, regionTestType, traitType) {
    m1 <- sum(mu.a * g.sum)
    var1 <- sum(mu.a * (1 - mu.a) * g.sum^2)
    var_b <- sum(Phi)
    q_adj <- m1 + (q.sum - m1) * sqrt(var1 / var_b)
    qinv <- -sign(q_adj - m1) * abs(q_adj - m1) + m1
    pval_noadj <- pchisq((q_adj - m1)^2 / var1, lower.tail = FALSE, df = 1, log.p = TRUE)
    is_spa_converge <- TRUE
    p_burden <- if (abs(q_adj - m1) / sqrt(var1) < 2) {
      pval_noadj
    } else {
      spa_pval(mu.a, g.sum, q_adj, qinv, pval_noadj, .Machine$double.eps^0.25, TRUE,
               traitType, is_spa_converge)
    }
    adj <- get("SPA_ER_kernel_related_Phiadj_fast_new", envir = ns)
    adj(p.new, Score, Phi, p_burden, regionTestType)
  }
  # Fix B: delete the burden-level rescale, i.e. the published SAIGE-GENE
  # adjustment (Supplementary Note 1.3.3): per-variant calibration only.
  adj_b <- edit_stmt(orig_adj, "r", "min(r, 1)", quote(r <- 1))
  list(skat103 = skat103, phi_a = phi_a, adj_b = adj_b,
       orig_skat = orig_skat, orig_phi = orig_phi, orig_adj = orig_adj)
}

# Install the fixes named in `spec` ("103", "A", "B" in any combination).
diag_install <- function(spec, fx) {
  if (grepl("103", spec)) assignInNamespace("get_SKAT_pvalue", fx$skat103, ns = "SAIGE")
  if (grepl("A", spec)) {
    assignInNamespace("get_newPhi_scaleFactor_traitType", fx$phi_a, ns = "SAIGE")
  }
  if (grepl("B", spec)) {
    assignInNamespace("SPA_ER_kernel_related_Phiadj_fast_new", fx$adj_b, ns = "SAIGE")
  }
  invisible(NULL)
}

# One SKAT-O cell recomputed offline from a phi_XX.rds dump under a variant
# spec, by the same namespace functions the pipeline runs.
diag_cell <- function(d, fx, spec, r_corr, region_test_type = "SKAT-O") {
  adj_fun <- if (grepl("B", spec)) fx$adj_b else fx$orig_adj
  phi_fun <- if (grepl("A", spec)) fx$phi_a else fx$orig_phi
  skat_fun <- if (grepl("103", spec)) fx$skat103 else fx$orig_skat
  assignInNamespace("SPA_ER_kernel_related_Phiadj_fast_new", adj_fun, ns = "SAIGE")
  on.exit(assignInNamespace("SPA_ER_kernel_related_Phiadj_fast_new", fx$orig_adj, ns = "SAIGE"))
  re <- phi_fun(d$q.sum, d$mu, d$g.sum, d$p.new, d$Score, d$Phi, region_test_type, "binary")
  out <- skat_fun(d$Score, re$val, r_corr, region_test_type)
  out$B_adj <- sum(re$val)
  out
}
