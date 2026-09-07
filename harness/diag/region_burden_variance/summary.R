# One compact table per (SKAT-O, BURDEN) pair from the analysis.rds files.
repo <- Sys.getenv("SAIGE_REPO", getwd())
out_dir <- Sys.getenv("DIAG_OUT", file.path(repo, "harness/work/diag-burden"))
options(width = 220, digits = 7)
labels <- c(
  "committed ratio file, default cutoffs (goldens P6/P5)",
  "every ratio 1.0, default cutoffs (the ledger's reference)",
  "sparse ratio 2.0, null 1.0: the clamp min(r, 1) is defused",
  "sparse ratio 0.5, null 1.0: the region output does not move",
  "committed ratios, boundary moved to MAC 200 (pseudo-variant on the sparse arm)",
  "sparse 2.0 AND boundary 200: both arms on the same GRM and ratio, no clamp",
  "null ratio 0.5, sparse 1.0: BURDEN halves, the region does not move"
)
for (k in 0:6) {
  d <- readRDS(file.path(out_dir, paste0("R", k), "analysis.rds"))
  cat(sprintf("\n===== pair R%d/B%d: %s =====\n", k, k, labels[k + 1]))
  print(data.frame(
    cell = d$cell, m = d$m, MAC = d$MAC, r = round(d$r, 5),
    `Badj/naive` = round(d$B_adj / d$var1_naive, 6),
    `Braw/naive` = round(d$B_raw / d$var1_naive, 6),
    `varB/naive` = round(d$var_B_implied / d$var1_naive, 6),
    `BETA_Rfix/BETA_B` = round(d$BETA_R_fix / d$BETA_B, 7),
    p_R = signif(d$p_R, 6), p_B = signif(d$p_B, 6), check.names = FALSE
  ), row.names = FALSE)
}
