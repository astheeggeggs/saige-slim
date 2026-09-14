#!/usr/bin/env python3
"""Replicate quantitative nulls for the pooled-row calibration (2026-09-14).

The p2 sweep (handoff-2026-09-14-flexrv-p2.md, 4a) put the pooled Cauchy
rows of FlexRV AND of the standard 6-cell grid at 0.09-0.095 on the single
`y_null` at alpha 0.05, outside the binomial interval of 210 genes, with the
tail clean. One draw cannot say whether 0.09 is the number or the draw. This
writes R fresh null phenotypes with the fixture's own recipe
(02_build_dataset.py, section 8): a polygenic component transmitted through
the pedigree (founders N(0, h2_poly); each child the midparent plus
N(0, h2_poly / 2)), the age/sex covariate effect at var_covar, and Gaussian
noise, so var = 1 and chr20 h2 = 0 exactly. Own RNG stream per replicate.

Writes pheno_pcs_null.txt beside the fixture: every column of pheno_pcs.txt
plus y_null_r1 .. y_null_rR. usage: make_null_series.py --dir <fixture> [-R 5]
"""
import argparse, json, os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True)
ap.add_argument("-R", type=int, default=5)
ap.add_argument("--seed", type=int, default=20260915)
a = ap.parse_args()
meta = json.load(open(os.path.join(a.dir, "sim_meta.json")))
h2_poly = meta["h2_poly"]; var_covar = 0.10           # 02_build_dataset.py --var-covar default
fam = [l.split() for l in open(os.path.join(a.dir, "chr20_sim.fam"))]
iid = [r[1] for r in fam]; pos = {x: i for i, x in enumerate(iid)}
N = len(fam)
ph = open(os.path.join(a.dir, "pheno_pcs.txt")).read().splitlines()
hdr = ph[0].split()
assert [l.split()[hdr.index("IID")] for l in ph[1:]] == iid, "pheno rows are not in fam order"
age = np.array([float(l.split()[hdr.index("age")]) for l in ph[1:]])
sex = np.array([float(l.split()[hdr.index("sex")]) for l in ph[1:]])
z_age = (age - age.mean()) / age.std()
cov_raw = 0.7 * z_age + 1.0 * (sex - sex.mean())
cov = cov_raw * np.sqrt(var_covar / cov_raw.var())
founder = np.array([r[2] == "0" and r[3] == "0" for r in fam])
child = np.where(~founder)[0]
fa = np.array([pos[fam[i][2]] for i in child]); mo = np.array([pos[fam[i][3]] for i in child])
assert founder[fa].all() and founder[mo].all(), "a parent is not a founder: the recipe assumes one generation"
cols = []
for r in range(1, a.R + 1):
    rng = np.random.default_rng(a.seed + r)
    poly = np.zeros(N)
    poly[founder] = rng.normal(0, np.sqrt(h2_poly), founder.sum())
    poly[child] = 0.5 * (poly[fa] + poly[mo]) + rng.normal(0, np.sqrt(h2_poly / 2), len(child))
    y = poly + cov + rng.normal(0, np.sqrt(1 - h2_poly - var_covar), N)
    cols.append(y)
    print(f"y_null_r{r}: var {y.var():.3f} var(poly) {poly.var():.3f} (target {h2_poly}) "
          f"cor(child, midparent) {np.corrcoef(poly[child], 0.5 * (poly[fa] + poly[mo]))[0, 1]:.3f}")
with open(os.path.join(a.dir, "pheno_pcs_null.txt"), "w") as fh:
    fh.write(" ".join(hdr + [f"y_null_r{r}" for r in range(1, a.R + 1)]) + "\n")
    for i, l in enumerate(ph[1:]):
        fh.write(l + " " + " ".join(f"{c[i]:.6f}" for c in cols) + "\n")
print(f"wrote pheno_pcs_null.txt with {a.R} replicate nulls ({founder.sum()} founders, {len(child)} children)")
