# fix_check.R: the fixed SKAT-O runs against BURDEN mode (plan step 1) and the
# no-op checks (plan step 4), from the out.txt files fixruns.sh wrote.
repo <- Sys.getenv("SAIGE_REPO", getwd())
out_dir <- Sys.getenv("DIAG_OUT", file.path(repo, "harness/work/diag-burden"))
options(width = 200, digits = 7)
read_cells <- function(r) {
  d <- read.table(file.path(out_dir, r, "out.txt"), header = TRUE, sep = "\t",
                  check.names = FALSE)
  d[d$Group != "Cauchy", ]
}
cell_names <- function(d) paste(d$Region, d$Group, d$max_MAF)
pairs <- list(
  c("R103", "B0", "#103 only, boundary 20.5 (today's code otherwise)"),
  c("RA0", "B0", "#103 + fix A, boundary 20.5"),
  c("RB0", "B0", "#103 + fix B, boundary 20.5"),
  c("RA4", "B4", "#103 + fix A, boundary 200 (both modes on the sparse arm)"),
  c("RB4", "B4", "#103 + fix B, boundary 200 (both modes on the sparse arm)")
)
for (p in pairs) {
  r <- read_cells(p[1])
  b <- read_cells(p[2])
  stopifnot(all(r$MAC == b$MAC))
  cat(sprintf("\n===== %s vs %s: %s =====\n", p[1], p[2], p[3]))
  print(data.frame(
    cell = cell_names(r), MAC = r$MAC,
    `BETA_R/BETA_B` = round(r$BETA_Burden / b$BETA_Burden, 5),
    `SE_R/SE_B` = round(r$SE_Burden / b$SE_Burden, 5),
    p_R = signif(r$Pvalue_Burden, 5), p_B = signif(b$Pvalue_Burden, 5),
    `p_R/p_B` = round(r$Pvalue_Burden / b$Pvalue_Burden, 4),
    p_SKATO = signif(r$Pvalue, 5), p_SKAT = signif(r$Pvalue_SKAT, 5),
    check.names = FALSE
  ), row.names = FALSE)
}
cols <- c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden", "SE_Burden")
rel <- function(a, b) {
  vapply(cols, function(cn) abs(b[[cn]] - a[[cn]]) / abs(a[[cn]]), numeric(nrow(a)))
}
cat("\n== fix A vs fix B, per cell max relative difference over the five columns ==\n")
for (pr in list(c("RA0", "RB0"), c("RA4", "RB4"))) {
  a <- read_cells(pr[1])
  b <- read_cells(pr[2])
  m <- apply(rel(a, b), 1, max)
  cat(sprintf("%s vs %s:\n", pr[1], pr[2]))
  print(data.frame(cell = cell_names(a), max_rel_diff = signif(m, 3)), row.names = FALSE)
}
cat("\n== no-op check: fix A (with #103) on the sparse-ratio-2.0 file, where the clamp was\n",
    "already inert, vs the unpatched run -- p-value columns only (#103 moves BETA/SE) ==\n")
a <- read_cells("R2")
b <- read_cells("R2A")
pcols <- c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT")
m <- apply(vapply(pcols, function(cn) abs(b[[cn]] - a[[cn]]) / abs(a[[cn]]), numeric(nrow(a))),
           1, max)
print(data.frame(cell = cell_names(a), max_rel_pdiff = signif(m, 3),
                 note = ifelse(m < 1e-6, "unchanged", "moved")), row.names = FALSE)
