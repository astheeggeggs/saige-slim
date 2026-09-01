#!/usr/bin/env bash
# LD-pruned GRM set, variance-ratio set, PCs, and the PGEN mirror.
set -euo pipefail
WT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$WT"
# $P2/$RS come from 00_env.sh, which addresses the pixi env ABSOLUTELY. A path
# relative to this checkout breaks in a git worktree, which has no .pixi of its
# own -- that is a silent "No such file or directory" partway through the run.
. "$(dirname "$0")/00_env.sh"
# Prefer an inherited D so run_all.sh exporting either variable is sufficient.
D="${SIMDIR:-${D:-$WT/harness/work/sim}}"
PFX=$D/chr20_sim

# --nonfounders everywhere: the 1000 offspring have PAT/MAT set, so plink2
# would otherwise compute allele frequencies from founders only, which would
# not reconcile with the counts SAIGE sees over all 10,000 samples.

# 1. LD-pruned set for the GRM. The upstream recipe wants imputed ARRAY data;
#    there is none here, so the stand-in is the common (MAF>1%) WGS variants.
#    The recipe's window is "50 5" -- 50 VARIANTS. On array data (~1 marker per
#    4 kb) that spans ~200 kb and removes most LD; on 2.27M WGS variants it
#    spans 1410 bp and removes almost none. Measured: the 50-variant window
#    leaves 14,047 markers whose effective independent count is only 1,732
#    (GRM off-diagonal sd 0.0240), while a 1 Mb PHYSICAL window leaves 3,055
#    markers with M_eff 2,702 (sd 0.0192) -- 2.7x fewer spurious related pairs
#    from 4.6x fewer markers. Use the physical window.
$P2 --bfile $PFX --nonfounders --maf 0.01 --indep-pairwise 1000kb 0.05 --out $D/grm --silent
$P2 --bfile $PFX --extract $D/grm.prune.in --make-bed --out $D/data_for_grm --silent

# 2. variance-ratio set: randomly chosen variants above and below MAC 20,
#    matching step 1's default MAC categories (10, 20.5] and >20.5.
$P2 --bfile $PFX --nonfounders --freq counts --out $D/allcounts --silent
python3 - "$D" <<'PY'
import random, sys
d = sys.argv[1]
random.seed(11)
lo, hi = [], []
with open(f"{d}/allcounts.acount") as fh:
    hdr = fh.readline().split()
    ia, io, ii = hdr.index("ALT_CTS"), hdr.index("OBS_CT"), hdr.index("ID")
    for line in fh:
        f = line.split()
        alt, obs = float(f[ia]), float(f[io])
        mac = min(alt, obs - alt)
        if 11 <= mac <= 20:
            lo.append(f[ii])
        elif mac > 20:
            hi.append(f[ii])
sl = random.sample(lo, min(2000, len(lo)))
sh = random.sample(hi, min(2000, len(hi)))
with open(f"{d}/vr_markers.txt", "w") as out:
    out.write("\n".join(sl + sh) + "\n")
print(f"variance-ratio markers: {len(sl)} with MAC in [11,20] (pool {len(lo):,}), "
      f"{len(sh)} with MAC>20 (pool {len(hi):,})")
PY
$P2 --bfile $PFX --extract $D/vr_markers.txt --make-bed --out $D/data_for_vr --silent

# 3. PCs from the LD-pruned set.
#    Exact --pca builds the dense N x N GRM: 8*N^2 bytes, which is 800 MB at
#    N=10,000 and 20 GB at N=50,000 -- plink2 aborts with "Out of memory.
#    Consider --pca approx instead. Failed allocation size: 20000000000", which
#    is exactly 8*50000^2. plink2's own guidance is that approx is preferred
#    above a few thousand samples, so switch on N rather than on failure.
NSAMP=$(wc -l < $PFX.fam)
PCA_MODE=""
if [ "$NSAMP" -gt 5000 ]; then
  PCA_MODE="approx"
  echo "PCA: N=$NSAMP > 5000, using --pca approx (exact needs $(( 8 * NSAMP / 1000 * NSAMP / 1000000 )) GB)"
fi
$P2 --bfile $D/data_for_grm --nonfounders --pca 10 $PCA_MODE --out $D/pca --silent

# 4. PGEN mirror (step-2 pgen reader).
#    psam-cols=sex forces a "#IID SEX" psam. SAIGE's PgenClass::readPsamFile
#    (src/PGEN.cpp:248) skips the '#' header and takes line_elements[0] as the
#    sample ID -- that is the FID under plink2's DEFAULT column set, so a plain
#    --make-pgen produces a psam step 2 rejects with "At least one subject
#    requested is not in psam file". The bundled fixture has FID == IID so no
#    golden covers this. Workaround, not a fix (bug-for-bug rule).
$P2 --bfile $PFX --make-pgen psam-cols=sex --out $D/chr20_sim_pgen --silent

wc -l < $D/grm.prune.in  | xargs echo "LD-pruned markers:"
wc -l < $D/data_for_vr.bim | xargs echo "variance-ratio markers:"
