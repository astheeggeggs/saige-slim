# Does the sparse GRM recover the planted pedigree, and is SAIGE's default
# relatedness cutoff of 0.05 defensible on a single-chromosome GRM?
#
# Truth: 500 families of 4 (FA, MO, CH1, CH2); each contributes 4
# parent-offspring pairs + 1 sib pair (kinship 0.25, GRM ~0.5); FA-MO are
# unrelated. The other 8000 samples are unrelated founders.
suppressMessages(library(Matrix))
d <- Sys.getenv("D", "harness/work/sim")
mtx <- commandArgs(TRUE)[1]
if (is.na(mtx)) mtx <- Sys.glob(file.path(d, "sparseGRM*randomMarkersUsed.sparseGRM.mtx"))[1]
ids <- readLines(paste0(mtx, ".sampleIDs.txt"))
G <- as(readMM(mtx), "TsparseMatrix")
fam <- read.table(file.path(d, "chr20_sim.fam"), stringsAsFactors = FALSE)
stopifnot(identical(ids, fam$V2))
pos <- setNames(seq_along(ids), ids)

tp <- list()
for (f in unique(fam$V1[fam$V1 != "UNREL"])) {
  m <- fam[fam$V1 == f, ]
  tp[[length(tp) + 1]] <- rbind(c(m$V2[1], m$V2[3]), c(m$V2[1], m$V2[4]),
                                c(m$V2[2], m$V2[3]), c(m$V2[2], m$V2[4]),
                                c(m$V2[3], m$V2[4]))
}
tp <- do.call(rbind, tp)
key <- function(i, j) paste(pmin(i, j), pmax(i, j))
truth <- key(pos[tp[, 1]], pos[tp[, 2]])
off <- G@i != G@j
k <- key(G@i[off] + 1, G@j[off] + 1)
v <- G@x[off]
istrue <- k %in% truth
np <- choose(length(ids), 2) - length(truth)

cat(sprintf("sparse GRM %s\n  %d x %d, %d off-diagonal pairs, %d true pairs planted\n",
            basename(mtx), nrow(G), ncol(G), length(v), length(truth)))
cat(sprintf("  value on TRUE pairs : median %.3f (IQR %.3f-%.3f)   [expected 0.5]\n",
            median(v[istrue]), quantile(v[istrue], .25), quantile(v[istrue], .75)))
cat(sprintf("  value on OTHER pairs: median %.3f, 99.9%% %.3f, max %.3f\n",
            median(v[!istrue]), quantile(v[!istrue], .999), max(v[!istrue])))
cat(sprintf("  diagonal: median %.3f, range %.3f-%.3f\n\n",
            median(G@x[!off]), min(G@x[!off]), max(G@x[!off])))
# "mass" = summed off-diagonal weight. What matters for the model is not the
# COUNT of spurious pairs but how much off-diagonal covariance they inject
# relative to the real pedigree, and how far that pushes rows away from
# diagonal dominance (max off-diagonal row sum vs a diagonal of ~1).
cat(sprintf("%-8s %12s %9s %7s %10s %10s %9s\n", "cutoff", "extra_prs", "true_prs",
            "sens%", "extra/true", "mass_ratio", "maxRowSum"))
n <- length(ids)
cand <- c(0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.25, 0.30, 0.40)
ok <- rep(FALSE, length(cand))
for (ci in seq_along(cand)) {
  cut <- cand[ci]
  s <- v >= cut
  ex <- sum(s & !istrue); tr <- sum(s & istrue)
  mx <- sum(v[s & !istrue]); mt <- sum(v[s & istrue])
  ii <- c(G@i[off][s], G@j[off][s]) + 1L
  w <- tapply(c(v[s], v[s]), factor(ii, levels = seq_len(n)), sum)
  w[is.na(w)] <- 0
  cat(sprintf("%-8.3f %12d %9d %7.1f %10.1f %10.2f %9.2f\n", cut, ex, tr,
              100 * tr / length(truth), ex / max(tr, 1), mx / max(mt, 1e-9), max(w)))
  ok[ci] <- (mx / max(mt, 1e-9) <= 0.5) && (tr / length(truth) >= 0.95)
}
# Selection rule, fixed in advance: the smallest cutoff at which the spurious
# off-diagonal MASS is at most half the real pedigree mass, while still
# retaining >=95% of the planted relative pairs. SAIGE's default of 0.05 is
# calibrated for a genome-wide GRM (off-diagonal SD ~0.005); on one chromosome
# the SD is ~0.025 so 0.05 is only ~2 SD out and cannot be used.
sel <- if (any(ok)) cand[which(ok)[1]] else max(cand)
cat(sprintf("\nCHOSEN_CUTOFF %s%s\n", sel,
            if (any(ok)) "" else "  (no candidate met the rule; using the largest)"))
