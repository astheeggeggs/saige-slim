# Numeric distance of each SKAT-O run's output from R0's, per column.
repo <- Sys.getenv("SAIGE_REPO", getwd())
out_dir <- Sys.getenv("DIAG_OUT", file.path(repo, "harness/work/diag-burden"))
read_cells <- function(r) {
  d <- read.table(file.path(out_dir, r, "out.txt"), header = TRUE, sep = "\t",
                  check.names = FALSE)
  d[d$Group != "Cauchy", ]
}
a <- read_cells("R0")
cols <- c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden", "SE_Burden")
for (r in c("R1", "R3", "R4", "R6", "R2", "R5")) {
  b <- read_cells(r)
  rel <- vapply(cols, function(cn) max(abs(b[[cn]] - a[[cn]]) / abs(a[[cn]])), numeric(1))
  cat(sprintf("%s vs R0: max rel diff  %s\n", r,
              paste(sprintf("%s=%.1e", cols, rel), collapse = "  ")))
}
