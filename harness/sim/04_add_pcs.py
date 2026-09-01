#!/usr/bin/env python
"""Merge plink2 PCs into the phenotype file (SAIGE wants one flat table)."""
import sys

d = sys.argv[1]
npc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
pcs = {}
with open(f"{d}/pca.eigenvec") as fh:
    hdr = fh.readline().split()
    iid = hdr.index("IID")
    for line in fh:
        f = line.split()
        pcs[f[iid]] = f[iid + 1:iid + 1 + npc]
with open(f"{d}/pheno.txt") as fh, open(f"{d}/pheno_pcs.txt", "w") as out:
    hdr = fh.readline().split()
    j = hdr.index("IID")
    out.write(" ".join(hdr + [f"PC{k + 1}" for k in range(npc)]) + "\n")
    n = 0
    for line in fh:
        f = line.split()
        out.write(" ".join(f + pcs[f[j]]) + "\n")
        n += 1
print(f"pheno_pcs.txt: {n} samples, {npc} PCs")
