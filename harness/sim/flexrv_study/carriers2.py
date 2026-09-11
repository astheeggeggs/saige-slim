#!/usr/bin/env python3
"""The collapse bound (docs/design/flexrv-sim-study.md, question 3): per gene,
the fraction of individuals carrying >= 2 qualifying ultra-rare variants
(MAC <= --collapse) in the cell, straight from the .bed. SAIGE collapses
ultra-rare variants with max(w*G), the paper with additive pseudo-variants;
the two differ only on these individuals, so this bounds the gap a priori.

usage: carriers2.py --dir <fixture> --groupfile <score-arm group file>
                    [--classes high,medium] [--maxmaf 0.001] [--collapse 10]
Reads MAF from variants.tsv.gz (the fixture's own in-sample AF), the bed as
PLINK 1 packed bytes, one row per marker. Prints per-gene counts and the
overall summary; writes <out> as TSV.
"""
import argparse, gzip, os, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True)
ap.add_argument("--groupfile", required=True)
ap.add_argument("--classes", default="high,medium")
ap.add_argument("--maxmaf", type=float, default=0.001)
ap.add_argument("--collapse", type=float, default=10)
ap.add_argument("--prefix", default="chr20_sim")
ap.add_argument("--out", default="")
a = ap.parse_args()
cls_ok = set(a.classes.split(","))

# marker index and MAF from the bim + variants table
bim = [l.split()[1] for l in open(os.path.join(a.dir, a.prefix + ".bim"))]
idx = {v: i for i, v in enumerate(bim)}
n = sum(1 for _ in open(os.path.join(a.dir, a.prefix + ".fam")))
maf = {}
with gzip.open(os.path.join(a.dir, "variants.tsv.gz"), "rt") as fh:
    h = fh.readline().rstrip("\n").split("\t"); ii, im, ia = h.index("ID"), h.index("MAF"), h.index("AC")
    for l in fh:
        f = l.rstrip("\n").split("\t")
        maf[f[ii]] = (float(f[im]), int(f[ia]))
nbytes = (n + 3) // 4
bed = open(os.path.join(a.dir, a.prefix + ".bed"), "rb")
assert bed.read(3) == b"\x6c\x1b\x01", "not a PLINK 1 SNP-major bed"
# packed 2-bit codes -> dosage of allele 1 (bim col 5, the ALT here): 00=2, 10=1, 11=0, 01=missing->0
lut = np.zeros((256, 4), np.int8)
for b in range(256):
    for s in range(4):
        c = (b >> (2 * s)) & 3
        lut[b, s] = {0: 2, 2: 1, 3: 0, 1: 0}[c]

lines = open(a.groupfile).read().splitlines()
rows = []
i = 0
while i < len(lines):
    v = lines[i].split(); an = lines[i + 1].split()
    assert v[1] == "var" and an[1] == "anno"
    gene = v[0]
    sel = [x for x, c in zip(v[2:], an[2:]) if c in cls_ok and x in maf
           and maf[x][0] <= a.maxmaf and min(maf[x][1], 2 * n - maf[x][1]) <= a.collapse
           and maf[x][1] > 0]
    i += 2
    while i < len(lines) and lines[i].split()[1] != "var":
        i += 1
    if not sel:
        rows.append((gene, 0, 0, 0, 0.0)); continue
    load = np.zeros(n, np.int32)
    for x in sel:
        bed.seek(3 + idx[x] * nbytes)
        raw = np.frombuffer(bed.read(nbytes), np.uint8)
        g = lut[raw].ravel()[:n]
        # alt-first bim: the minor allele may be allele 2; count the MINOR allele
        if maf[x][1] * 2 > 2 * n:   # AC counts the ALT; if ALT is the major allele flip
            g = 2 - g
        load += g
    n1 = int((load >= 1).sum()); n2 = int((load >= 2).sum())
    rows.append((gene, len(sel), n1, n2, n2 / n1 if n1 else 0.0))

tot_ur = sum(r[1] for r in rows); tot1 = sum(r[2] for r in rows); tot2 = sum(r[3] for r in rows)
print(f"genes {len(rows)}, ultra-rare markers in the cell {tot_ur}, carriers {tot1}, "
      f">=2-carriers {tot2}  ({tot2 / max(tot1, 1) * 100:.2f}% of carriers; "
      f"{tot2 / n / len(rows) * 100:.3f}% of individuals per gene on average)")
fr = np.array([r[4] for r in rows if r[2] > 0])
if len(fr):
    print("per-gene fraction of carriers with >= 2: quantiles 0/50/90/100 =",
          np.round(np.quantile(fr, [0, .5, .9, 1]), 4))
if a.out:
    with open(a.out, "w") as fh:
        fh.write("gene\tn_ur\tcarriers\tcarriers_ge2\tfrac_ge2\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.4f}\n")
