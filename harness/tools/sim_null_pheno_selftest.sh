#!/bin/bash
# Self-test for sim_null_pheno.R. Controls before conclusions, each on the bundled fixture
# (1,000 SIMULATED samples: counts are written plainly) or on a six-sample FABRICATED model:
#   1. reproducible under a seed, different under another;
#   2. positive control, quantitative: over six seeds the residual variance reproduces
#      tau1 + tau2 x mean(diag K) and the GRM pair covariance reproduces tau2 x k_ij -- on the
#      eigen path AND on the sparse-Cholesky path forced by --maxDenseComponent=2;
#   3. negative control: --tau2=0 draws no random effect, residual variance = tau1;
#   4. binary: the simulated prevalence is the fixed part's;
#   5. a planted INDEFINITE family block is clipped on the eigen path (count, min eigenvalue
#      equal to R's own) and jittered on the Cholesky path, both reported; a PD block is neither;
#   6. a phenotype row outside the model gets NA, and --mask masks a small people count;
#   7. step 1 refits on the simulated file with the golden command and converges.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "${HERE}/../.." && pwd)
source "${REPO}/harness/env.sh" >/dev/null 2>&1 || true
R="${RSCRIPT:-${REPO}/.pixi/envs/default/bin/Rscript}"
T=$(mktemp -d); trap 'rm -rf "${T}"' EXIT
G="${REPO}/extdata/output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
P="${REPO}/extdata/input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt"
MQ="${REPO}/harness/goldens/step1/quant_sparseGRM.rda"; MB="${REPO}/extdata/output/example_binary_sparseGRM.rda"
run() { # tag model phenoCol extra...
  local tag=$1 model=$2 col=$3; shift 3
  "${R}" "${HERE}/sim_null_pheno.R" --modelFile="${model}" --sparseGRMFile="${G}" --sparseGRMSampleIDFile="${G}.sampleIDs.txt" \
    --phenoFile="${P}" --phenoCol="${col}" --sampleIDCol=IID --outPhenoFile="${T}/${tag}.txt" --reportFile="${T}/${tag}.tsv" --mask=0 "$@" > "${T}/${tag}.log" 2>&1 \
    || { echo "FAIL ${tag}: $(tail -3 "${T}/${tag}.log")"; exit 1; }
}
val() { awk -F'\t' -v k="$2" '$1==k{print $2}' "$1"; }
fail=0; nchk=0
ok() { nchk=$((nchk+1)); echo "ok   $1"; }
bad() { fail=1; echo "FAIL $1"; }

# 1. seeds
run s1a "${MQ}" y_quantitative --seed=1 & run s1b "${MQ}" y_quantitative --seed=1 & run s2 "${MQ}" y_quantitative --seed=2 & wait
cmp -s "${T}/s1a.txt" "${T}/s1b.txt" && ok "seed 1 reproduces byte for byte" || bad "seed 1 not reproducible"
cmp -s "${T}/s1a.txt" "${T}/s2.txt" && bad "seed 2 equals seed 1" || ok "seed 2 differs"

# 2. positive control, both paths, six seeds each (in parallel)
for s in 1 2 3 4 5 6; do run e${s} "${MQ}" y_quantitative --seed=${s} & run c${s} "${MQ}" y_quantitative --seed=${s} --maxDenseComponent=2 & done; wait
for path in e c; do
  python3 - "${T}" "${path}" <<'PY' && ok "quantitative moments on path ${path}: residual variance and pair covariance reproduce tau" || bad "quantitative moments on path ${path}"
import sys
T, path = sys.argv[1], sys.argv[2]
def rep(f): return dict(l.rstrip("\n").split("\t") for l in open(f) if "\t" in l and not l.startswith("field"))
rs = [rep(f"{T}/{path}{s}.tsv") for s in range(1, 7)]
tau1, tau2, md = float(rs[0]["tau1"]), float(rs[0]["tau2"]), float(rs[0]["grm_mean_diag"])
v = sum(float(r["simulated_var_resid"]) for r in rs) / 6; exp = tau1 + tau2 * md
pc = sum(float(r["pair_cov_observed"]) / float(r["pair_cov_expected"]) for r in rs) / 6
print(f"     path {path}: mean var_resid {v:.4f} vs {exp:.4f}; pair cov ratio {pc:.3f}")
assert abs(v / exp - 1) < 0.10, f"var_resid {v} vs {exp}"
assert 0.75 < pc < 1.25, f"pair cov ratio {pc}"
if path == "c": assert all(r["components_jittered"] == "0" for r in rs), "a PD block needed jitter"
PY
done

# 3. negative control
run z "${MQ}" y_quantitative --seed=3 --tau2=0
python3 -c "import sys; assert float(sys.argv[1]) == 0" "$(val "${T}/z.tsv" var_b)" && ok "tau2=0 draws no random effect" || bad "tau2=0 var_b $(val "${T}/z.tsv" var_b)"
python3 -c "import sys; t1=float(sys.argv[1]); v=float(sys.argv[2]); assert abs(v/t1-1)<0.15, (v,t1)" "$(val "${T}/z.tsv" tau1)" "$(val "${T}/z.tsv" simulated_var_resid)" && ok "tau2=0 residual variance = tau1" || bad "tau2=0 residual variance"

# 4. binary
for s in 1 2 3; do run b${s} "${MB}" y_binary --seed=${s} & done; wait
python3 - "${T}" <<'PY' && ok "binary prevalence is the fixed part's (three seeds)" || bad "binary prevalence"
import sys; T = sys.argv[1]
def rep(f): return dict(l.rstrip("\n").split("\t") for l in open(f) if "\t" in l and not l.startswith("field"))
rs = [rep(f"{T}/b{s}.tsv") for s in (1, 2, 3)]
p = sum(float(r["simulated_prevalence"]) for r in rs) / 3; e = float(rs[0]["fixed_part_prevalence"])
print(f"     prevalence {p:.4f} vs fixed part {e:.4f}, observed {rs[0]['observed_prevalence']}")
assert abs(p / e - 1) < 0.2, (p, e)
PY
grep -q "^rows_with_value" "${T}/b1.tsv" && [[ $(awk 'NR>1{print $2}' "${T}/b1.txt" | sort -u | tr '\n' ' ') == "0 1 " ]] && ok "binary column holds 0/1" || bad "binary column"

# 5. fabricated six-sample model: block A (3) planted indefinite, block B (2) PD, one singleton
"${R}" - "${T}" <<'RS'
T <- commandArgs(TRUE)[1]
ids <- paste0("s", 1:6)
modglmm <- list(theta = c(0.5, 1.0), coefficients = 0.1, offset = c(0.2, -0.1, 0, 0.3, -0.2, 0.1), isCovariateOffset = TRUE,
                sampleID = ids, traitType = "quantitative", y = rnorm(6))
save(modglmm, file = file.path(T, "fake.rda"))
# block A: s1 related to s2 and s3 at 0.9, s2-s3 absent (a thinned entry): eigenvalues
# 1 +- 0.9 sqrt(2), min -0.2728 -- indefinite. The pd twin has 0.4: min 1 - 0.4 sqrt(2) = 0.434.
for (f in c("ind", "pd")) {
  k <- if (f == "ind") .9 else .4
  A <- matrix(c(1, k, k, k, 1, 0, k, 0, 1), 3); cat(f, "block min eig:", min(eigen(A, symmetric = TRUE)$values), "\n")
  e <- rbind(cbind(1:6, 1:6, 1), c(1, 2, k), c(1, 3, k), c(4, 5, .5))
  writeLines(c("%%MatrixMarket matrix coordinate real symmetric", paste(6, 6, nrow(e)),
               apply(e, 1, function(r) paste(r[2], r[1], r[3]))), file.path(T, paste0(f, ".mtx")))
}
writeLines(ids, file.path(T, "fake.ids"))
write.table(data.frame(IID = c(ids, "extra"), y = 0, x1 = 1), file.path(T, "fake.pheno"), quote = FALSE, row.names = FALSE)
RS
for f in ind pd; do for cap in 2000 2; do
  "${R}" "${HERE}/sim_null_pheno.R" --modelFile="${T}/fake.rda" --sparseGRMFile="${T}/${f}.mtx" --sparseGRMSampleIDFile="${T}/fake.ids" \
    --phenoFile="${T}/fake.pheno" --phenoCol=y --sampleIDCol=IID --outPhenoFile="${T}/${f}${cap}.out" --reportFile="${T}/${f}${cap}.tsv" --mask=0 --maxDenseComponent=${cap} --seed=5 > "${T}/${f}${cap}.log" 2>&1 \
    || { echo "FAIL fake ${f} cap ${cap}: $(tail -2 "${T}/${f}${cap}.log")"; fail=1; }
done; done
[[ $(val "${T}/ind2000.tsv" components) == 3 && $(val "${T}/ind2000.tsv" components_eigen_clipped) == 1 ]] && ok "planted indefinite block: 3 components, 1 clipped on the eigen path" || bad "indefinite eigen: $(cat "${T}/ind2000.tsv" | tr '\n' ';')"
python3 -c "import sys; a=float(sys.argv[1]); assert abs(a-(-0.2727922))<1e-5, a" "$(val "${T}/ind2000.tsv" min_eigenvalue)" && ok "min eigenvalue -0.27279 reported (1 - 0.9 sqrt 2)" || bad "min eigenvalue $(val "${T}/ind2000.tsv" min_eigenvalue)"
[[ $(val "${T}/ind2.tsv" components_jittered) == 1 ]] && python3 -c "import sys; assert 0.27 < float(sys.argv[1]) <= 0.5" "$(val "${T}/ind2.tsv" max_jitter)" && ok "planted indefinite block jittered on the Cholesky path (jitter $(val "${T}/ind2.tsv" max_jitter), above the 0.27 deficit)" || bad "indefinite cholesky: $(cat "${T}/ind2.tsv" | tr '\n' ';')"
[[ $(val "${T}/pd2000.tsv" components_eigen_clipped) == 0 && $(val "${T}/pd2.tsv" components_jittered) == 0 ]] && ok "PD blocks neither clipped nor jittered" || bad "pd blocks"

# 6. NA outside the model; masking
[[ $(awk 'NR>1 && $1=="extra"{print $2}' "${T}/pd2000.out") == "NA" ]] && ok "a row outside the model is NA" || bad "extra row: $(grep extra "${T}/pd2000.out")"
"${R}" "${HERE}/sim_null_pheno.R" --modelFile="${T}/fake.rda" --sparseGRMFile="${T}/pd.mtx" --sparseGRMSampleIDFile="${T}/fake.ids" \
  --phenoFile="${T}/fake.pheno" --phenoCol=y --sampleIDCol=IID --outPhenoFile="${T}/m.out" --reportFile="${T}/m.tsv" --mask=20 > "${T}/m.log" 2>&1
[[ $(val "${T}/m.tsv" component_max_size) == "<20" && $(val "${T}/m.tsv" components_singleton) == "<20" ]] && ok "--mask=20 masks people counts below 20" || bad "mask: $(grep -E "max_size|singleton" "${T}/m.tsv" | tr '\n' ';')"

# 7. step 1 on the simulated quantitative file (the golden's command)
( cd "${REPO}/extdata" && "${R}" step1_fitNULLGLMM.R \
  --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
  --sparseGRMFile="${G}" --sparseGRMSampleIDFile="${G}.sampleIDs.txt" --useSparseGRMtoFitNULL=TRUE \
  --phenoFile="${T}/e1.txt" --phenoCol=y_quantitative --covarColList=x1,x2,a9,a10 --qCovarColList=a9,a10 \
  --sampleIDColinphenoFile=IID --traitType=quantitative --invNormalize=TRUE --IsOverwriteVarianceRatioFile=TRUE \
  --isCateVarianceRatio=TRUE --outputPrefix="${T}/step1" ) > "${T}/step1.log" 2>&1 || { bad "step 1 died: $(tail -3 "${T}/step1.log")"; }
"${R}" -e "load('${T}/step1.rda'); stopifnot(isTRUE(modglmm\$converged)); cat('     step-1 theta', modglmm\$theta, '\n')" > "${T}/s1chk.log" 2>&1 && { cat "${T}/s1chk.log"; ok "step 1 refits on the simulated file and converges"; } || bad "step 1 model: $(cat "${T}/s1chk.log")"

(( fail )) && { echo "SELFTEST FAILED"; exit 1; }
echo "selftest OK: ${nchk} checks"
