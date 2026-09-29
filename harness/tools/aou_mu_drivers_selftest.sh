#!/bin/bash
# Self-test of aou_mu_drivers.sh on simulated models, where the driver is known.
#   A  the bundled pinned fixture: the separating covariate `sep` carries the tail
#   B  a synthetic cohort whose risk falls steeply with age (a quadratic in
#      log-odds, the shape of a condition of reproductive age in an older
#      cohort): `age` and `age2` carry the tail, a bystander covariate does not
#   C  a synthetic cohort with a `sex` column that is CONSTANT: step 1 drops a
#      collinear column before the fit, so the model has no such column
#   D  a model with no tail: the summary is withheld, the covariates are listed
#   E  masking: a tail of fewer than MASK samples is written <MASK and withheld
#   F  a design saved without column names (found by timing the tool on a
#      full-size synthetic model): decomposed under numbered columns
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
RS="${RSCRIPT:-${REPO}/.pixi/envs/default/bin/Rscript}"
W="$(mktemp -d)"; trap 'rm -rf "${W}"' EXIT
mkdir -p "${W}/m"
G="${REPO}/harness/goldens/step1"
cp "${G}/binary_sparseGRM_pinned.rda" "${W}/m/A_pinned.rda"
cp "${G}/binary_sparseGRM_a.rda" "${W}/m/D_notail.rda"
"${RS}" -e '
  w <- commandArgs(TRUE)[1]; set.seed(175); n <- 6000
  age <- runif(n, 20, 85); pc1 <- rnorm(n); bystander <- rbinom(n, 1, 0.4)
  eta <- -2.2 - 0.004 * (age - 30)^2 + 0.1 * pc1                    # 1e-6 and below past ~85
  mu <- plogis(eta); y <- rbinom(n, 1, mu)
  qrx <- function(X) { Q <- qr.Q(qr(X)) * sqrt(nrow(X)); colnames(Q) <- colnames(X); Q }
  X <- cbind(minus1 = -1, age = age, age2 = age^2, bystander1 = bystander, PC1 = pc1)
  modglmm <- list(fitted.values = mu, y = y, X = qrx(X), traitType = "binary")
  save(modglmm, file = file.path(w, "B_age.rda"))
  # C: the covariates as they are AFTER step 1 dropped a constant sex column
  modglmm <- list(fitted.values = mu, y = y, X = qrx(X[, c("minus1", "age", "age2", "PC1")]), traitType = "binary")
  save(modglmm, file = file.path(w, "C_nosex.rda"))
  # E: a tail of 7 samples
  mu2 <- plogis(-2 + 0.1 * pc1); mu2[1:7] <- 1e-9; y2 <- rbinom(n, 1, mu2); y2[1:7] <- 0
  modglmm <- list(fitted.values = mu2, y = y2, X = qrx(X), traitType = "binary")
  save(modglmm, file = file.path(w, "E_small.rda"))
  # F: the design saved WITHOUT column names
  Xu <- qrx(X); colnames(Xu) <- NULL
  modglmm <- list(fitted.values = mu, y = y, X = Xu, traitType = "binary")
  save(modglmm, file = file.path(w, "F_unnamed.rda"))
  cat("B: tail at 1e-6:", sum(ifelse(y == 1, 1 - mu, mu) <= 1e-6), "\n")' "${W}/m"
fails=0; checks=0
check(){ checks=$((checks+1)); if [[ "$2" == "$3" ]]; then echo "  ok   $1"; else echo "  FAIL $1: got '$2', want '$3'"; fails=$((fails+1)); fi; }
cd "${W}"
MODELS="${W}/m" OUT_DIR="${W}/o" TAG=t MASK=20 RSCRIPT="${RS}" bash "${REPO}/harness/tools/aou_mu_drivers.sh" > "${W}/log" 2>&1 || { cat "${W}/log"; exit 1; }
gunzip -c "${W}/o/mu_drivers.t.tsv.gz" > "${W}/t.tsv"
top(){ # top <model>: the covariate with the largest |tail_minus_rest|
  awk -F'\t' -v M="$1" 'NR>1 && $1==M && $8!="withheld"{v=$8+0; if(v<0)v=-v; if(v>b){b=v; n=$5}} END{print n}' "${W}/t.tsv"; }
sharegap(){ # sharegap <model> <regex>: percent of the summed gap carried by the covariates matching
  awk -F'\t' -v M="$1" -v R="$2" 'NR>1 && $1==M && $8!="withheld"{t+=$8; if($5 ~ R) s+=$8} END{printf "%d", (t==0)?0:100*s/t+0.5}' "${W}/t.tsv"; }
check "six binary models read" "$(grep -c '^===== ' "${W}/log")" "6"
check "F: an unnamed design is decomposed under numbered columns" "$(top F_unnamed.rda | sed 's/[23]$//')" "column"
check "A: the separating covariate carries the tail" "$(top A_pinned.rda)" "sep1"
a=$(sharegap A_pinned.rda '^sep1$'); check "A: it carries more than 90 % of the gap" "$(( a > 90 ))" "1"
check "B: age or its square is the largest driver" "$(top B_age.rda | sed 's/2$//')" "age"
b=$(sharegap B_age.rda '^age2?$'); check "B: age and age2 together carry more than 90 % of the gap" "$(( b > 90 ))" "1"
c=$(sharegap B_age.rda '^bystander1$'); check "B: the bystander carries under 5 % of the gap" "$(( c < 5 && c > -5 ))" "1"
check "B: the summary names no sex column" "$(grep -A4 '^===== B_age.rda' "${W}/log" | grep -c "columns whose name contains 'sex': none")" "1"
check "C: a dropped constant column is absent from the covariates" "$(awk -F'\t' '$1=="C_nosex.rda" && $5 ~ /sex/' "${W}/t.tsv" | grep -c . || true)" "0"
check "A: a column named for the separating covariate is listed" "$(grep -A3 '^===== A_pinned.rda' "${W}/log" | grep -c 'covariates in the model: x1 x2 a91 a92 a101 a102 sep1')" "1"
check "D: no tail, statistics withheld" "$(grep -A6 '^===== D_notail.rda' "${W}/log" | grep -c 'tail statistics withheld')" "1"
check "E: a tail of 7 is written <20" "$(awk -F'\t' '$1=="E_small.rda"{print $4}' "${W}/t.tsv" | sort -u)" "<20"
check "E: its rows are withheld" "$(awk -F'\t' '$1=="E_small.rda"{print $8}' "${W}/t.tsv" | sort -u)" "withheld"
check "E: the log does not print the 7" "$(grep -A2 '^===== E_small.rda' "${W}/log" | grep -c ': 7 *$' || true)" "0"
# the decomposition is exact: per model the covariates' gaps and the remainder sum to the gap the log states
for mdl in A_pinned.rda B_age.rda; do
  s=$(awk -F'\t' -v M="${mdl}" 'NR>1 && $1==M{t+=$8} END{printf "%.2f", t}' "${W}/t.tsv")
  g=$(grep -A5 "^===== ${mdl}" "${W}/log" | grep "linear predictor sits" | awk '{printf "%.2f", $6}')
  check "${mdl}: the parts sum to the stated gap (${g})" "$(awk -v a="${s}" -v b="${g}" 'BEGIN{d=a-b; if(d<0)d=-d; print (d <= 0.02*(b<0?-b:b)+0.02) ? "yes" : "no"}')" "yes"
done
echo "${checks} checks, ${fails} failed"
exit $(( fails > 0 ))
