#!/usr/bin/env python
"""Write the sparse GRM implied by the pedigree, exactly.

Stage 3. `createSparseGRM` estimates relatedness from genotypes; over one region
it cannot do so well, because the number of effectively independent common
markers is bounded by the region's genetic length rather than by how many
markers it contains (25 Mb of chr20 is ~41 cM). The pedigree is known here, so
the answer can be written down instead of estimated, and step 0 becomes a
measured CONTROL -- score its output against this -- rather than an unchecked
dependency of steps 1 and 2.

Output matches what `Matrix::writeMM` emits for SAIGE's own sparse GRM, because
`setSparseSigma_new` (R/SAIGE_SPATest_Region.R:26) reads it with
`Matrix:::readMM`: a `coordinate real symmetric` file holding the LOWER triangle,
1-indexed, whose indices address the companion `.sampleIDs.txt` line by line.

Values are on SAIGE's GRM scale, which is twice the kinship: 1.0 on the diagonal
and 0.5 for a first-degree pair. Founders are independent coalescent draws and
the two parents of a family are unrelated, so the within-family first-degree
pairs are the complete set of nonzero off-diagonals -- there is nothing else to
include, and that is a property of how 02_build_dataset.py builds the pedigree,
not an approximation.
"""
import argparse
import os

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True, help="dataset directory (02's --outdir)")
ap.add_argument("--prefix", default="chr20_sim")
ap.add_argument("--out", default=None,
                help="output .mtx (default <dir>/pedigree_truth.sparseGRM.mtx)")
a = ap.parse_args()

out = a.out or os.path.join(a.dir, "pedigree_truth.sparseGRM.mtx")
fam = os.path.join(a.dir, f"{a.prefix}.fam")
ped = os.path.join(a.dir, f"{a.prefix}.pedigree_truth.txt")
for f in (fam, ped):
    if not os.path.exists(f):
        raise SystemExit(f"missing {f} -- run 02_build_dataset.py first")

# .fam order IS the sample order of the .bed, so the GRM must use it too.
iids = [l.split()[1] for l in open(fam)]
index = {s: i + 1 for i, s in enumerate(iids)}      # 1-indexed, as Matrix Market
assert len(index) == len(iids), "duplicate IIDs in .fam"
N = len(iids)

entries = [(i, i, 1.0) for i in range(1, N + 1)]    # diagonal
n_pair = 0
with open(ped) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ci = {c: k for k, c in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[ci["relation"]] == "self":
            continue
        x, y = index[f[ci["IID1"]]], index[f[ci["IID2"]]]
        assert x != y, f
        # lower triangle: row > col
        entries.append((max(x, y), min(x, y), float(f[ci["grm"]])))
        n_pair += 1

entries.sort(key=lambda e: (e[0], e[1]))            # row-major, as SAIGE emits
assert len(entries) == N + n_pair
assert all(r >= c for r, c, _ in entries)

with open(out, "w") as fh:
    fh.write("%%MatrixMarket matrix coordinate real symmetric\n")
    fh.write(f"{N} {N} {len(entries)}\n")
    for r, c, v in entries:
        fh.write(f"{r} {c} {v}\n")
with open(out + ".sampleIDs.txt", "w") as fh:
    fh.write("\n".join(iids) + "\n")

print(f"wrote {out}")
print(f"  {N} x {N}, {len(entries)} stored entries "
      f"({N} diagonal + {n_pair} first-degree pairs)")
print(f"  density {len(entries) / (N * (N + 1) / 2) * 100:.4f}% of the lower triangle")
print(f"wrote {out}.sampleIDs.txt ({N} IDs, .fam order)")
