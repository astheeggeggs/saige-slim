#!/usr/bin/env python3
"""Build step-2 marker fixtures with a CONTROLLED fraction of ER-eligible markers.

    pixi run python3 harness/make_rare_fixture.py <out-dir> <f> [<f> ...]

Why this exists
---------------
The retirement decision between the two step-2 parallelism mechanisms
(docs/design/step2-threading-verdict.md) turns on one number that the bundled
fixture cannot vary: the fraction f of markers that are ER-eligible, i.e.
binary with MAC <= max_MAC_use_ER.  Those markers are MAIN-THREAD ONLY under
in-process threading -- Rf_lchoose aborts off the main thread -- so they are
decoded on a worker, discarded, and re-run serially.  Fork has no such limit.
On the bundled fixture f = 0.041, which is far too small to discriminate: both
mechanisms measured identical there (p = 0.14 / 0.53 / 0.49).

The obvious cheap proxy -- raise --max_MAC_use_ER to make more markers eligible
-- DOES NOT WORK, and it is worth saying why rather than rediscovering it.  ER
cost climbs steeply with MAC, and above ~21 carriers the exact test switches to
the resampling arm entirely.  Raising the cutoff therefore changes the WORKLOAD
(each deferred marker becomes far more expensive) instead of changing only f.
The comparison would no longer be like for like.

So f is varied here the honest way: by changing which markers are in the file,
holding the per-marker cost of an ER marker at its realistic value.

What it does NOT do
-------------------
It does not simulate genotypes.  Every marker written is a byte-for-byte copy of
a real marker's genotype block from the bundled 22-chromosome fixture, drawn
from one of two pools:

    rare pool    real markers with 1 <= MAC <= 4   (ER-eligible; 5,292 of them)
    common pool  real markers with MAC > 4         (not eligible)

For a target f the output has round(f*M) markers cycled from the rare pool and
the rest cycled from the common pool, with M and N matching the original so the
existing step-1 model applies unchanged and the serial baseline stays
comparable across f.

Two choices that matter, and would bias the result if made differently:

* The two kinds are INTERLEAVED UNIFORMLY, not concatenated.  parallelFor hands
  out contiguous ranges, so a file with all its rare markers in one block would
  dump the entire serial burden on one worker and measure scheduling rather
  than the mechanism.  Uniform interleaving gives every range the same f, which
  is also what a real genome looks like.
* Positions are RENUMBERED consecutively.  SPAGMMATtest sorts markers by
  CHROM/POS, so duplicated markers carrying their original positions would be
  reordered relative to file order and the interleaving would be destroyed.

Markers are reused (a pool of 5,292 cannot fill 128,868 slots otherwise).  That
is sound for a timing fixture -- decode and test cost per marker do not depend
on whether another marker shares its genotypes -- but it means these files are
for BENCHMARKING ONLY.  They are not goldens and nothing about their p-values
is meaningful.
"""
import os, sys

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "extdata", "input",
                   "nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr")
MAC_ER = 4          # matches the --max_MAC_use_ER default the benchmark uses

# PLINK1 .bed 2-bit codes, low bits first within each byte:
#   00 hom A1, 01 missing, 10 het, 11 hom A2
_A2 = []; _MISS = []
for b in range(256):
    a2 = miss = 0
    for i in range(4):
        c = (b >> (2 * i)) & 3
        if c == 3: a2 += 2
        elif c == 2: a2 += 1
        elif c == 1: miss += 1
    _A2.append(a2); _MISS.append(miss)


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__.strip().splitlines()[2].strip())
    outdir, targets = sys.argv[1], [float(x) for x in sys.argv[2:]]
    os.makedirs(outdir, exist_ok=True)

    fam = open(SRC + ".fam", "rb").read()
    n = fam.count(b"\n")
    bim = open(SRC + ".bim").read().splitlines()
    m = len(bim)
    per = (n + 3) // 4
    bed = open(SRC + ".bed", "rb").read()
    assert bed[:3] == b"\x6c\x1b\x01", "not a SNP-major PLINK1 .bed"
    assert len(bed) == 3 + m * per, (len(bed), 3 + m * per)
    print(f"source: {m} markers x {n} samples, {per} bytes/marker")

    rare, common = [], []
    for j in range(m):
        blk = bed[3 + j * per: 3 + (j + 1) * per]
        a2 = miss = 0
        for b in blk:
            a2 += _A2[b]; miss += _MISS[b]
        nn = n - miss
        mac = min(a2, 2 * nn - a2)
        if mac == 0:
            continue                      # monomorphic: --minMAC=1 would drop it
        (rare if mac <= MAC_ER else common).append(j)
    print(f"pools: rare (MAC<={MAC_ER}) {len(rare)}, common {len(common)}")

    for f in targets:
        n_rare = int(round(f * m))
        tag = f"rare{int(round(f * 100)):03d}"
        out = os.path.join(outdir, tag)
        # Uniform interleave: marker i is rare iff floor((i+1)*f) > floor(i*f).
        pick, ri, ci = [], 0, 0
        for i in range(m):
            if int((i + 1) * n_rare / m) > int(i * n_rare / m):
                pick.append(rare[ri % len(rare)]); ri += 1
            else:
                pick.append(common[ci % len(common)]); ci += 1
        with open(out + ".bed", "wb") as fh:
            fh.write(b"\x6c\x1b\x01")
            for j in pick:
                fh.write(bed[3 + j * per: 3 + (j + 1) * per])
        with open(out + ".bim", "w") as fh:
            for i, j in enumerate(pick):
                c = bim[j].split("\t")
                # chrom 1, consecutive positions, unique IDs: keeps file order ==
                # sorted order, so the interleaving survives SPAGMMATtest's sort.
                fh.write("\t".join(["1", f"m{i+1}", "0", str(i + 1), c[4], c[5]]) + "\n")
        with open(out + ".fam", "wb") as fh:
            fh.write(fam)
        print(f"  {tag}: {len(pick)} markers, {ri} rare ({100.0*ri/m:.1f}%) -> {out}.bed")


if __name__ == "__main__":
    main()
