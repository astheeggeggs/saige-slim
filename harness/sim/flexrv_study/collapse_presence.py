#!/usr/bin/env python3
"""Which (region, weight set) rows exist at collapse 0 but not at collapse 10, or
vice versa, over every FlexRV pair in a results directory. usage:
collapse_presence.py RESULTS_DIR. A row present at c0 and absent at c10 is a set
whose entire support is ultra-rare and whose weighted pseudo-marker was zeroed
by --dosage_zerod_cutoff (src/UTIL.cpp imputeGenoAndFlip: MAC in WEIGHTED units
<= --dosage_zerod_MAC_cutoff, entries <= the cutoff cleaned); at c0 the same
markers stay individual 0/1/2 genotypes and survive."""
import sys, os, glob, collections
R = sys.argv[1]
def keys(fn):
    out = {}
    with open(fn) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        iw = hdr.index("Weight") if "Weight" in hdr else None
        im = hdr.index("MAC")
        for l in f:
            c = l.rstrip("\n").split("\t")
            if c[1] == "Cauchy" or iw is None: continue
            out[(c[0], c[iw])] = c[im]
    return out
tot0, tot10 = collections.Counter(), collections.Counter()
for f0 in sorted(glob.glob(f"{R}/*_c0_r1.txt")):
    f10 = f0.replace("_c0_r1", "_c10_r1")
    if not os.path.exists(f10) or os.path.getsize(f0) == 0 or os.path.getsize(f10) == 0: continue
    k0, k10 = keys(f0), keys(f10)
    if not k0: continue
    only0 = set(k0) - set(k10); only10 = set(k10) - set(k0)
    fam0 = collections.Counter(w.split("_", 1)[1] if "_" in w else w for _, w in only0)
    print(f"{os.path.basename(f0).replace('_c0_r1.txt', ''):40s} sets c0 {len(k0):6d} c10 {len(k10):6d}  only-c0 {len(only0):4d}  only-c10 {len(only10):3d}  "
          f"MAC of only-c0 rows: {collections.Counter(k0[k] for k in only0).most_common(4)}")
    for w, n in fam0.items(): tot0[w] += n
    for _, w in only10: tot10[w] += 1
print("\nonly-at-c0 by set (all files):", sorted(tot0.items(), key=lambda x: -x[1])[:20])
print("only-at-c10 by set (all files):", sorted(tot10.items(), key=lambda x: -x[1])[:10])
