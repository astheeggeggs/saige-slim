#!/bin/bash
# Where do a null model's fitted values sit? One row per binary model: how many
# samples have a FLIP PROBABILITY (mu for a control, 1 - mu for a case -- the
# probability of the outcome the sample does not have) in each decade, and the
# mass below each candidate threshold of the pinned-sample guard (LEDGER #174).
#
# WHY. The guard's threshold was set at the fit gate's 1e-8 on the assumption
# that a separated group is fitted below it. The bundled fixture says it is
# not: R's glm stops on the relative change of the deviance, which leaves a
# separated group near 1e-8 x deviance / (1.26 x group size) -- 1.7e-9..1.2e-7
# for 60 samples of 1,000, with 13 of the 60 below 1e-8 and 47 above. The step-2
# logs report only the count below 1e-8 (aou_pinned_scan.sh), so they cannot say
# how many samples sit just above it. This reads the models.
#
# Runs INSIDE a dsub job on the SAIGE image (it needs R, nothing else). The
# models are localized by dsub (--input MODELS=gs://.../nullglmm_v4/<anc>_*.rda),
# so nothing leaves the workspace; the table holds counts per decade with
# every count from 1 to MASK-1 written as "<MASK" and the masses that go with
# a masked count written as "masked". Safe to read and to keep. Env:
#   MODELS   a directory searched recursively for *.rda, OR a file pattern such
#            as /path/eas_*.rda -- which is what dsub sets the variable to for
#            `--input MODELS=gs://.../nullglmm_v4/eas_*.rda`, the form to use
#            when one directory holds every ancestry's models. Per-chromosome
#            LOCO files `_chr<j>.rda` are skipped: the model's own file is read.
#   OUT_DIR  local dir dsub delocalizes (--output-recursive)
#   TAG      label for the output file, default models
#   MASK     default 20
#   RSCRIPT  default Rscript
# Output: <OUT_DIR>/mu_decades.<TAG>.tsv.gz and the same table on stdout.
set -euo pipefail
: "${MODELS:?}"; : "${OUT_DIR:?}"; TAG="${TAG:-models}"; MASK="${MASK:-20}"; RSCRIPT="${RSCRIPT:-Rscript}"
mkdir -p "${OUT_DIR}"
if [[ -d ${MODELS} ]]; then root="${MODELS}"; pat='*.rda'; else root="$(dirname "${MODELS}")"; pat="$(basename "${MODELS}")"; fi
[[ -d ${root} ]] || { echo "ERROR: no such directory: ${root} (MODELS=${MODELS})"; exit 1; }
find "${root}" -type f -name "${pat}" -name '*.rda' ! -name '*_chr[0-9]*.rda' | sort > models.list
nfound=$(awk 'END{print NR}' models.list)
echo "models found: ${nfound} (directory ${root}, pattern ${pat})"
# nothing to read is an error, not an empty table: a wrong path must not look like a clean result
(( nfound > 0 )) || { echo "ERROR: no model file matched"; exit 1; }
cat > mu_decades.R <<'RS'
args <- commandArgs(trailingOnly = TRUE)
files <- readLines(args[1]); out <- args[2]; mask <- as.integer(args[3])
edges <- c(-Inf, -15, -12, -10, -9, -8, -7, -6, -5, -4, -3, 0)
labs <- c("le1e-15", "1e-15_1e-12", "1e-12_1e-10", "1e-10_1e-9", "1e-9_1e-8", "1e-8_1e-7",
          "1e-7_1e-6", "1e-6_1e-5", "1e-5_1e-4", "1e-4_1e-3", "gt1e-3")
thr <- c(1e-8, 1e-7, 1e-6, 1e-5)
show <- function(k) ifelse(k > 0 & k < mask, paste0("<", mask), as.character(k))
rows <- list()
for (f in files) {
  e <- new.env()
  ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(err) FALSE)
  if (!ok || is.null(e$modglmm)) { cat("skipped (no modglmm):", basename(f), "\n"); next }
  m <- e$modglmm
  if (!identical(as.character(m$traitType), "binary")) next
  mu <- as.numeric(m$fitted.values); y <- as.numeric(m$y)
  flip <- ifelse(y == 1, 1 - mu, mu)
  # a fitted value at the boundary whose outcome CONTRADICTS it is not pinned: counted apart
  contra <- sum((mu <= 1e-6 & y == 1) | (mu >= 1 - 1e-6 & y == 0))
  k <- as.integer(table(cut(log10(pmax(flip, 1e-300)), edges, labels = labs)))
  below <- vapply(thr, function(t) sum(flip <= t), numeric(1))
  mass <- vapply(thr, function(t) sum(flip[flip <= t]), numeric(1))
  masked <- below > 0 & below < mask
  # the widest run of empty decades above the lowest occupied one: a separated
  # group shows as a gap between it and the rest
  occ <- which(k > 0)
  gap <- if (length(occ) > 1) max(diff(occ)) - 1L else 0L
  r <- c(model = basename(f), N = length(mu), cases = show(sum(y == 1)),
         stats::setNames(show(k), labs),
         stats::setNames(show(below), paste0("n_le_", format(thr, scientific = TRUE))),
         stats::setNames(ifelse(masked, "masked", formatC(mass, format = "e", digits = 2)),
                         paste0("mass_le_", format(thr, scientific = TRUE))),
         empty_decades_in_widest_gap = gap, contradicting = show(contra))
  rows[[length(rows) + 1]] <- r
}
if (length(rows) == 0) { cat("no binary model read\n"); quit(status = 1) }
tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
utils::write.table(tab, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("binary models read:", nrow(tab), "\n")
n8 <- sum(tab[["n_le_1e-08"]] != "0"); n6 <- sum(tab[["n_le_1e-06"]] != "0"); n5 <- sum(tab[["n_le_1e-05"]] != "0")
cat("models with any sample at or below 1e-8:", n8, "; 1e-6:", n6, "; 1e-5:", n5, "\n")
cat("models whose samples at or below 1e-6 are not all at or below 1e-8:",
    sum(tab[["n_le_1e-06"]] != tab[["n_le_1e-08"]]), "\n")
RS
out="${OUT_DIR}/mu_decades.${TAG}.tsv"
"${RSCRIPT}" mu_decades.R models.list "${out}" "${MASK}"
cat "${out}"
gzip -f "${out}"
echo done
