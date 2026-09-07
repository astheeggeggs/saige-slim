# Instrumented step-2 driver for the region burden-variance diagnosis
# (docs/state/diag-2026-09-07-region-burden-variance.md). Optionally installs
# candidate fixes (SAIGE_FIX, see fixes.R), then wraps three R-side region
# functions in the SAIGE namespace so every call's inputs and outputs are
# dumped as .rds files under $SAIGE_DUMP, then sources the unmodified CLI
# script. Nothing in the package is touched; without SAIGE_FIX the output
# files stay byte-identical to an uninstrumented run (goldens P5/P6/P36/P37).
suppressPackageStartupMessages(library(SAIGE))
ns <- asNamespace("SAIGE")
repo <- Sys.getenv("SAIGE_REPO")
fix <- Sys.getenv("SAIGE_FIX", "")
if (nzchar(fix)) {
  source(file.path(repo, "harness", "diag", "region_burden_variance", "fixes.R"))
  diag_install(fix, diag_fixes(ns))
}
dump <- Sys.getenv("SAIGE_DUMP", "")
if (nzchar(dump)) {
  orig_phi <- get("get_newPhi_scaleFactor_traitType", envir = ns)
  orig_adj <- get("SPA_ER_kernel_related_Phiadj_fast_new", envir = ns)
  orig_skat <- get("get_SKAT_pvalue", envir = ns)
  counter <- new.env()
  counter$phi <- 0L
  counter$adj <- 0L
  counter$skat <- 0L
  wrap_phi <- function(q.sum, mu.a, g.sum, p.new, Score, Phi, regionTestType, traitType) {
    out <- orig_phi(q.sum, mu.a, g.sum, p.new, Score, Phi, regionTestType, traitType)
    counter$phi <- counter$phi + 1L
    saveRDS(list(q.sum = q.sum, mu = mu.a, g.sum = g.sum, p.new = p.new,
                 Score = Score, Phi = Phi, out = out),
            file.path(dump, sprintf("phi_%02d.rds", counter$phi)))
    out
  }
  wrap_adj <- function(p.new, Score, Phi, p_value_burden, regionTestType) {
    out <- orig_adj(p.new, Score, Phi, p_value_burden, regionTestType)
    counter$adj <- counter$adj + 1L
    saveRDS(list(p.new = p.new, Score = Score, Phi = Phi,
                 p_value_burden = p_value_burden, out = out),
            file.path(dump, sprintf("adj_%02d.rds", counter$adj)))
    out
  }
  wrap_skat <- function(Score, Phi, r.corr, regionTestType) {
    out <- orig_skat(Score, Phi, r.corr, regionTestType)
    counter$skat <- counter$skat + 1L
    saveRDS(list(Score = Score, Phi = Phi, r.corr = r.corr, out = out),
            file.path(dump, sprintf("skat_%02d.rds", counter$skat)))
    out
  }
  assignInNamespace("get_newPhi_scaleFactor_traitType", wrap_phi, ns = "SAIGE")
  assignInNamespace("SPA_ER_kernel_related_Phiadj_fast_new", wrap_adj, ns = "SAIGE")
  assignInNamespace("get_SKAT_pvalue", wrap_skat, ns = "SAIGE")
}
source(file.path(repo, "extdata", "step2_SPAtests.R"))
