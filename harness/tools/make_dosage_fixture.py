#!/usr/bin/env python3
"""Build a DOSAGE pgen over the step-2 sample set, from the hard-call .bed.

`--dosage_zerod_cutoff` / `--dosage_zerod_MAC_cutoff` only do their real job on
fractional dosages: `t_GVec.clean(cutoff)` (src/UTIL.cpp:108-112) zeroes every
entry below the cutoff whenever the marker's MAC is at or below
`--dosage_zerod_MAC_cutoff`. On hard calls the only sub-threshold values are the
ones `--impute_method=mean` writes (LEDGER #109); on dosages it reaches the DATA.
The pgen dosage decode path itself (`PgrGet1D` + `Dosage16ToDoubles`,
src/PGEN.cpp:485-489) is likewise touched by no golden.

`extdata/input/dosage_20markers_100samples.*` already exists and carries real
dosages, but it has 100 samples and is wired only to test-step1-pgen.R's
hard-call rounding check. The bundled step-1 models carry 1,000 samples, so a
step-2 dosage golden needs a dosage pgen over THAT sample set -- hence this,
which keeps genotype_100markers' 10,000 samples and 100 markers so the result is
directly comparable to goldens P1 and P42.

Construction, deterministic and minimal:

  * every sample keeps its exact hard call (0, 1 or 2) as a dosage, EXCEPT
  * twelve hom-ref samples per marker, chosen by the stride `(j*89 + t*37) mod
    N_MODEL` and skipped if not hom-ref, which take a fixed ladder straddling
    the 0.2 default cutoff:
        0.02 0.05 0.10 0.15 0.19 0.199 | 0.21 0.25 0.30 0.50 0.80 1.30
    six strictly below and six at or above.

Only hom-ref samples are overwritten, so a marker's MAC can only RISE, and it
rises by at most sum(ladder) = 4.069. That is what keeps the low-MAC markers
below `--dosage_zerod_MAC_cutoff` (10) where the clean actually fires, while the
common markers stay above it -- the contrast the golden pair is for. The twelve
are taken from the MODEL's first N_MODEL samples for the reason recorded in
make_missing_fixture.py: step 2 analyses that subset, so anything placed outside
it is invisible.

NO MISSING CALLS: this fixture varies one thing. Missingness is
genotype_100markers_missing (goldens P42-P48).
"""
import sys

LADDER = [0.02, 0.05, 0.10, 0.15, 0.19, 0.199, 0.21, 0.25, 0.30, 0.50, 0.80, 1.30]
N_MODEL = 1000
MAGIC = b"\x6c\x1b\x01"
# PLINK-1 .bed codes, alt-first (A1 == ALT): 0 = hom A1, 1 = missing, 2 = het,
# 3 = hom A2. So the A1 (ALT) dosage is 2, ., 1, 0 respectively.
A1_DOSAGE = {0: 2.0, 2: 1.0, 3: 0.0}


def main(src, dst):
    n = sum(1 for _ in open(src + ".fam"))
    bim = [l.split() for l in open(src + ".bim")]
    raw = open(src + ".bed", "rb").read()
    assert raw[:3] == MAGIC
    body, m, per = raw[3:], len(bim), (n + 3) // 4

    rows, stats = [], []
    for j in range(m):
        off = j * per
        g = [(body[off + (i >> 2)] >> (2 * (i & 3))) & 3 for i in range(n)]
        assert 1 not in g, "source .bed has missing calls; this fixture is dosage-only"
        d = [A1_DOSAGE[c] for c in g]
        hard_mac_sub = sum(d[:N_MODEL])
        placed, t = [], 0
        for step in range(N_MODEL):
            if len(placed) == len(LADDER):
                break
            i = (j * 89 + step * 37) % N_MODEL
            if g[i] == 3 and i not in placed:        # hom-ref only: MAC can only rise
                d[i] = LADDER[len(placed)]
                placed.append(i)
        assert len(placed) == len(LADDER), f"marker {bim[j][1]}: only {len(placed)} hom-ref slots"
        dos_mac_sub = sum(d[:N_MODEL])
        stats.append((bim[j][1], hard_mac_sub, dos_mac_sub))
        # A1 is the ALT allele (bim column 5), A2 the REF (column 6)
        rows.append([bim[j][1], bim[j][4], bim[j][5]] +
                    [("%g" % v) for v in d])

    with open(dst + ".dosage", "w") as fh:
        for r in rows:
            fh.write(" ".join(r) + "\n")
    with open(dst + ".map", "w") as fh:
        for b in bim:
            fh.write(f"{b[0]}\t{b[1]}\t{b[2]}\t{b[3]}\n")
    for ext in (".fam",):
        open(dst + ext, "w").write(open(src + ext).read())

    fires = [s for s in stats if 1 <= s[2] <= 10]
    above = [s for s in stats if s[2] > 10]
    print(f"wrote {dst}.dosage  n={n} m={m}  (analysed subset = first {N_MODEL})")
    print(f"  ladder sum per marker                       : {sum(LADDER):.3f}")
    print(f"  markers with dosage MAC <= 10 (clean FIRES) : {len(fires)}")
    print(f"  markers with dosage MAC >  10 (clean inert) : {len(above)}")
    print(f"  example fires: {[s[0] for s in fires[:8]]}")
    assert len(fires) >= 5 and len(above) >= 5, "need both sides of the MAC gate"


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
