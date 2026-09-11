#!/bin/bash
# FlexRV simulation study driver (docs/design/flexrv-sim-study.md).
# One step-2 run per (trait, score arm, method, collapse); every row records
# commit, dirty flag, SAIGE.so md5, load and the fixture's seed. Frozen step-1
# models from the fixture directory; the pedigree-truth sparse GRM.
#
# usage: run.sh <trait> <arm> <method> <collapse> [rep]
#   trait     y_null | y_null_bin | y_quant | y_binary   (step1_<trait>.rda in $D)
#   arm       informative | realistic | uninformative | synonymous | plain
#   method    flex192 | flexscore | flexmaf | flex1 | unw | beta | skato | grid
#             (grid = the standard practice: SKAT-O over high,high:medium x
#              maxMAF 1e-4,1e-3,1e-2, six cells, the Cauchy row is the result)
#   collapse  0 | 10   (--MACCutoff_to_CollapseUltraRare)
# env: D (fixture dir), G (group-file dir from make_scores.py), OUT (results dir),
#      NT (threads, default 1), CELL (annotation group, default high:medium),
#      MAXMAF (default 0.001), MINLOAD_REFUSE (default 2.5)
set -euo pipefail
REPO=/Users/dpalmer/Repositories/SAIGE
D=${D:-/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k}
G=${G:-$REPO/harness/work/flexrv_study/scores}
OUT=${OUT:-$REPO/harness/work/flexrv_study/results}
NT=${NT:-1}; CELL=${CELL:-high:medium}; MAXMAF=${MAXMAF:-0.001}
trait=$1; arm=$2; method=$3; collapse=$4; rep=${5:-1}
if [ "$method" = grid ]; then CELL=high,high:medium; MAXMAF=0.0001,0.001,0.01; fi
mkdir -p "$OUT"
load=$(uptime | sed 's/.*load averages: //' | cut -d' ' -f1)
if awk -v l="$load" -v m="${MINLOAD_REFUSE:-2.5}" 'BEGIN{exit !(l > m)}'; then
  echo "REFUSE: load $load > ${MINLOAD_REFUSE:-2.5}" >&2; exit 3
fi
cd $REPO; source harness/env.sh
RS=$REPO/.pixi/envs/default/lib/R/bin/Rscript
gf=$G/groupFile_$arm.txt; [ "$arm" = plain ] && gf=$D/groupFile.txt
[ -f "$gf" ] || { echo "no group file $gf" >&2; exit 2; }
[ "$arm" = synonymous ] && CELL=low   # the low-only mask; COMMON reads CELL below
COMMON="--bedFile=$D/chr20_sim.bed --bimFile=$D/chr20_sim.bim --famFile=$D/chr20_sim.fam \
 --chrom=20 --minMAF=0 --minMAC=0.5 --LOCO=FALSE --is_fastTest=FALSE --nThreads=$NT \
 --GMMATmodelFile=$D/step1_$trait.rda --varianceRatioFile=$D/step1_$trait.varianceRatio.txt \
 --sparseGRMFile=$D/pedigree_truth.sparseGRM.mtx \
 --sparseGRMSampleIDFile=$D/pedigree_truth.sparseGRM.mtx.sampleIDs.txt \
 --groupFile=$gf --annotation_in_groupTest=$CELL --MACCutoff_to_CollapseUltraRare=$collapse"
FLEX="--flexRV_score=SIM --flexRV_maxMAF=$MAXMAF --flexRV_lofAnno=high"
# The synonymous arm is the low-only mask: no LoF label exists there, so the
# LoF score transform (all zeros, refused by the CLI) is dropped and the arm is
# 15 x 12 = 180 sets; a method's own --flexRV_scoreTransforms (later) wins.
if [ "$arm" = synonymous ]; then
  FLEX="$FLEX --flexRV_scoreTransforms=unw;cbrt:0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9;sbeta:1.5,2,3,5,10,20"
fi
case $method in
  flex192)   ARGS="$COMMON $FLEX --r.corr=1" ;;
  flexscore) ARGS="$COMMON $FLEX --r.corr=1 --flexRV_MAFTransforms=unw" ;;
  flexmaf)   ARGS="$COMMON $FLEX --r.corr=1 --flexRV_scoreTransforms=unw" ;;
  flex1)     ARGS="$COMMON $FLEX --r.corr=1 --flexRV_scoreTransforms=unw --flexRV_MAFTransforms=unw" ;;
  unw)       ARGS="$COMMON --maxMAF_in_groupTest=$MAXMAF --is_no_weight_in_groupTest=TRUE --r.corr=1" ;;
  beta)      ARGS="$COMMON --maxMAF_in_groupTest=$MAXMAF --r.corr=1" ;;
  skato)     ARGS="$COMMON --maxMAF_in_groupTest=$MAXMAF" ;;
  grid)      ARGS="$COMMON --maxMAF_in_groupTest=$MAXMAF" ;;
  *) echo "unknown method $method" >&2; exit 2 ;;
esac
tag=${trait}_${arm}_${method}_c${collapse}_r${rep}
out=$OUT/$tag.txt; log=$OUT/$tag.log
{
  echo "tag=$tag start=$(date '+%F %T') load=$load nthreads=$NT cell=$CELL maxmaf=$MAXMAF"
  echo "head=$(git rev-parse --short HEAD) dirty=$(git status --porcelain | grep -c .) so_md5=$(md5 -q .pixi/envs/default/lib/R/library/SAIGE/libs/SAIGE.so) seed=$(python3 -c "import json;print(json.load(open('$D/sim_meta.json'))['seed'])")"
  cd extdata
  /usr/bin/time -l $RS step2_SPAtests.R $ARGS --SAIGEOutputFile=$out
  echo "end=$(date '+%F %T') exit=$?"
} > $log 2>&1
n=$(grep -c . $out 2>/dev/null || echo 0)
echo "$tag: $n rows, $(grep -E ' real ' $log | awk '{print $3 " s user"}')"
