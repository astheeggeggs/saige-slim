#!/usr/bin/env bash
# Upstream SAIGE (main, 365deeff) single-variant step 2 on the SAME fixture,
# same model file, same GRM, single-threaded -- so the only variable is the
# implementation. main exposes PGEN as --pgenPrefix (one prefix); the strip
# branch takes --pgenFile/--pvarFile/--psamFile. Same capability, different CLI.
set -uo pipefail
M=/Users/dpalmer/Repositories/SAIGE-main
D=/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k
OUT=$M/bench; mkdir -p "$OUT"
export R_LIBS_USER=$M/Rlib
export RCPP_PARALLEL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
RS=$M/.pixi/envs/default/lib/R/bin/Rscript
GRM=$D/pedigree_truth.sparseGRM.mtx
REPS="${REPS:-3}"

common=(
  --chrom=20 --minMAF=0 --minMAC=20
  --GMMATmodelFile=$D/step1_y_quant.rda
  --varianceRatioFile=$D/step1_y_quant.varianceRatio.txt
  --sparseGRMFile=$GRM --sparseGRMSampleIDFile=$GRM.sampleIDs.txt
  --LOCO=FALSE --is_Firth_beta=TRUE --pCutoffforFirth=0.1 --is_fastTest=TRUE
  --nThreads=1
)
warm(){ cat "$D/chr20_sim.bed" >/dev/null 2>&1; cat "$D/chr20_sim_pgen.pgen" >/dev/null 2>&1; }

run_arm(){
  local arm=$1
  local rep=$2
  local log="$OUT/${arm}_${rep}.log"
  local res="$OUT/${arm}_${rep}.txt"
  local -a io
  if [ "$arm" = bed ]; then
    io=(--bedFile=$D/chr20_sim.bed --bimFile=$D/chr20_sim.bim --famFile=$D/chr20_sim.fam)
  elif [ "$arm" = pgen ]; then
    # Upstream REFUSES pgen without this: "genotype is in pgen, please set
    # AlleleOrder=ref-first". The strip branch's reader has no such requirement
    # and matches its own bed output byte for byte with no flag. Allele order
    # flips which allele is counted, so it changes BETA's sign and AF -- hence
    # the matched bedref arm below, to compare results at equal orientation.
    io=(--pgenPrefix=$D/chr20_sim_pgen --AlleleOrder=ref-first)
  else
    io=(--bedFile=$D/chr20_sim.bed --bimFile=$D/chr20_sim.bim --famFile=$D/chr20_sim.fam
        --AlleleOrder=ref-first)
  fi
  warm
  cd "$M"
  /usr/bin/time -p $RS extdata/step2_SPAtests.R "${io[@]}" "${common[@]}" \
      --SAIGEOutputFile="$res" > "$log" 2>&1
  local rc=$?
  echo -e "$arm\t$rep\t$rc\t$(awk '/^real/{print $2}' "$log"|tail -1)\t$(awk '/^user/{print $2}' "$log"|tail -1)"
}
echo -e "arm\trep\trc\treal_s\tuser_s"
for r in $(seq 1 $REPS); do run_arm bed $r; run_arm pgen $r; done
