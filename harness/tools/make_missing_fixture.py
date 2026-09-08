#!/usr/bin/env python3
"""Plant a deterministic missing-call pattern in a copy of a PLINK .bed.

Every bundled .bed has EXACTLY ZERO missing calls -- 7 files, ~150 million
calls, verified by decoding them -- so `--maxMissing`, `--impute_method`,
`--dosage_zerod_cutoff` and `--dosage_zerod_MAC_cutoff` have no golden and no
test. This writes the fixture that reaches them.

THE DENOMINATOR IS THE MODEL'S SAMPLE COUNT, NOT THE .fam's. Step 2 subsets the
genotype file to the null model's samples, and `missingRate = numMissing / N`
is computed over THAT subset (src/PLINK.cpp:311, src/PGEN.cpp:584). The bundled
step-1 models carry 1,000 samples and they are fam rows 1..1000 exactly, while
genotype_100markers.fam has 10,000 rows. A first version of this script sized
the ladder against 10,000 and planted by a stride over all of them; SAIGE then
reported rates of 0.032/0.096/0.161 where 0.05/0.10/0.15 were intended, the
boundary marker fell on the wrong side of `--maxMissing`, and the two impute
methods differed on 2 markers instead of 6. Size against N_MODEL.

    j % 10 :  0     1      2      3     4     5     6      7      8     9
    calls  :  0     1      5      10    50    100   150    151    300   600
    rate   :  0.0   0.001  0.005  0.01  0.05  0.10  0.15   0.151  0.30  0.60
    at the 0.15 default:  <-------- passes -------->  <----- dropped ----->

`j % 10 == 6` sits EXACTLY on the cutoff and must PASS, because the gate is a
strict `missingRate > g_cfg.missingRate_cutoff` (src/Main.cpp:468 and :531).
150/1000 and the parsed literal 0.15 are the same double, so this pins the
strict inequality rather than measuring floating-point luck.

Two properties the ladder must have -- checked below, not assumed, because a
golden that cannot differ from its twin gates nothing. Both are computed over
the MODEL's subset for the same reason as above:

  * `--impute_method` best_guess and minor differ ONLY where imputeG =
    round(2*AF) is 1 rather than 0, i.e. AF >= 0.25 (src/UTIL.cpp:83-97). Ten
    markers qualify over the first 1,000 samples (over all 10,000 it is nine,
    and a different nine -- another reason the subset matters).
  * `--dosage_zerod_cutoff` fires only when MAC <= `--dosage_zerod_MAC_cutoff`
    (default 10), where it zeroes what `mean` imputation just wrote
    (src/UTIL.cpp:108-112).

Which SAMPLES go missing, inside the subset, is the stride `(j*97 + i*7) mod
N_MODEL`; 7 is coprime with 1,000 so the k indices are distinct with no dedup
pass, and the start moves with j so missingness is not nested across markers.

OUT-OF-SUBSET CONTROL. Every marker also gets `OUT_OF_SUBSET` missing calls in
rows N_MODEL.. , which no bundled model analyses. They must NEVER reach the
output: the reported missing rate has to stay exactly k/N_MODEL. That is a live
check on sample subsetting in both readers -- the defect class of LEDGER #52 --
and it is free. It does mean the FILE's overall missing rate is not the reported
one; that is deliberate.
"""
import sys

LADDER_CALLS = [0, 1, 5, 10, 50, 100, 150, 151, 300, 600]   # per 1,000 samples
N_MODEL = 1000
OUT_OF_SUBSET = 500
MAGIC = b"\x6c\x1b\x01"


def calls(body, off, n):
    return [(body[off + (i >> 2)] >> (2 * (i & 3))) & 3 for i in range(n)]


def set_missing(body, off, s):
    byte, sh = off + (s >> 2), 2 * (s & 3)
    body[byte] = (body[byte] & ~(3 << sh)) | (1 << sh)      # 0b01 == missing


def main(src, dst):
    n = sum(1 for _ in open(src + ".fam"))
    ids = [l.split()[1] for l in open(src + ".bim")]
    raw = open(src + ".bed", "rb").read()
    assert raw[:3] == MAGIC, "not a variant-major PLINK-1 .bed"
    body, m, per = bytearray(raw[3:]), len(ids), (n + 3) // 4
    assert len(body) == m * per and n > N_MODEL

    stats = []
    for j in range(m):
        off, k = j * per, LADDER_CALLS[j % 10]
        sub = calls(body, off, N_MODEL)
        assert 1 not in sub, f"marker {ids[j]} already has missing calls in the subset"
        a1 = 2 * sub.count(0) + sub.count(2)
        af, mac = a1 / (2 * N_MODEL), min(a1, 2 * N_MODEL - a1)
        for i in range(k):
            set_missing(body, off, (j * 97 + i * 7) % N_MODEL)
        for i in range(OUT_OF_SUBSET):                       # never analysed
            set_missing(body, off, N_MODEL + ((j * 31 + i * 11) % (n - N_MODEL)))
        stats.append((j, ids[j], af, mac, k))

    passing = [s for s in stats if s[4] > 0 and s[4] <= 150]
    hi = [s for s in passing if s[2] >= 0.25]
    lo = [s for s in passing if 1 <= s[3] <= 10]
    at_cut = [s for s in stats if s[4] == 150]
    assert len(hi) >= 3, f"only {len(hi)} AF>=0.25 markers below the cutoff"
    assert len(lo) >= 3, f"only {len(lo)} MAC<=10 markers below the cutoff"
    assert len([s for s in stats if s[4] == 0]) >= 5
    assert at_cut, "no marker sits exactly on the --maxMissing default"

    with open(dst + ".bed", "wb") as fh:
        fh.write(MAGIC); fh.write(bytes(body))
    for ext in (".bim", ".fam"):
        open(dst + ext, "w").write(open(src + ext).read())

    print(f"wrote {dst}.bed  n={n} m={m}  (analysed subset = first {N_MODEL})")
    print(f"  untouched in-subset (byte-identity control): "
          f"{len([s for s in stats if s[4] == 0])}")
    print(f"  above --maxMissing=0.15 (dropped)          : "
          f"{len([s for s in stats if s[4] > 150])}")
    print(f"  EXACTLY at 0.15 (must pass, strict >)      : {[s[1] for s in at_cut]}")
    print(f"  AF>=0.25 and passing (best_guess != minor) : {len(hi)}  {[s[1] for s in hi]}")
    print(f"  MAC<=10 and passing (dosage_zerod fires)   : {len(lo)}  {[s[1] for s in lo]}")
    print(f"  out-of-subset decoys per marker            : {OUT_OF_SUBSET}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
