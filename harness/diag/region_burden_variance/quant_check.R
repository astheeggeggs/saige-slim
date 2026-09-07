# Quantitative negative control: no R-side Phi adjustment, no sparse GRM.
repo <- Sys.getenv("SAIGE_REPO", getwd())
out_dir <- Sys.getenv("DIAG_OUT", file.path(repo, "harness/work/diag-burden"))
options(width = 200, digits = 7)
read_cells <- function(r) {
  d <- read.table(file.path(out_dir, r, "out.txt"), header = TRUE, sep = "\t",
                  check.names = FALSE)
  d[d$Group != "Cauchy", ]
}
rt <- read_cells("RQ")
bt <- read_cells("BQ")
rows <- lapply(seq_len(nrow(rt)), function(k) {
  s <- readRDS(file.path(out_dir, "RQ", sprintf("skat_%02d.rds", k)))
  phi <- as.matrix(s$Phi)
  sc <- s$Score
  data.frame(
    cell = paste(rt$Region[k], rt$Group[k], rt$max_MAF[k]), m = length(sc), MAC = rt$MAC[k],
    `BETA_R==S/trPhi` = s$out$BETA_Burden / (sum(sc) / sum(diag(phi))),
    `BETA_Rfix/BETA_B` = (sum(sc) / sum(phi)) / bt$BETA_Burden[k],
    p_R = rt$Pvalue_Burden[k], p_B = bt$Pvalue_Burden[k], check.names = FALSE
  )
})
cat("\n== quantitative model (goldens P37 / P36): no adjustment, no sparse GRM ==\n")
print(do.call(rbind, rows), row.names = FALSE)
