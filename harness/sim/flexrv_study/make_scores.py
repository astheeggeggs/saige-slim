#!/usr/bin/env python3
"""FlexRV simulation study: the three score arms as SAIGE group files.

Reads a fixture's groupFile.txt (var + anno lines) and variants.tsv.gz (which
carries a_lat/e_lat/true_class since 2026-09-10) and writes one group file per
arm with a `score:SIM` line after each gene's anno line:

  informative   high -> 1; medium/low -> within-gene percentile rank of a_lat
  realistic     the same on a_lat + N(0, s^2), s = --noise (default 1.5), so the
                score is a noisy readout of the latent the effect class follows
  uninformative the informative scores permuted within gene (cor with the
                true class ~0): the negative control of the power question
  synonymous    the `low` class only, with scores resampled with replacement
                from the gene's `medium` scores (the paper's synonymous check;
                run it with --annotation_in_groupTest=low)

Own seed; the fixture is not touched. Prints, per arm, cor(score, true_class)
over the medium/low variants and over all coding variants, which is the
knob's calibration. Usage: make_scores.py --dir <fixture> --out <dir> [--noise s]
"""
import argparse, gzip, os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--noise", type=float, default=1.5)
ap.add_argument("--seed", type=int, default=20260910)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
rng = np.random.default_rng(a.seed)

alat, tcls = {}, {}
with gzip.open(os.path.join(a.dir, "variants.tsv.gz"), "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ia, ic, ii = hdr.index("a_lat"), hdr.index("true_class"), hdr.index("ID")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[ia] != ".":
            alat[f[ii]] = float(f[ia])
            tcls[f[ii]] = int(f[ic])

lines = open(os.path.join(a.dir, "groupFile.txt")).read().splitlines()
assert len(lines) % 2 == 0
arms = ["informative", "realistic", "uninformative", "synonymous"]
out = {k: [] for k in arms}
stat = {k: [] for k in arms}          # (score, true_class, is_high)


def pct(x):
    """within-gene percentile rank in (0, 1]; ties by average rank"""
    r = np.argsort(np.argsort(x, kind="stable"), kind="stable") + 1.0
    return r / len(x)


for i in range(0, len(lines), 2):
    v = lines[i].split(); an = lines[i + 1].split()
    assert v[1] == "var" and an[1] == "anno" and v[0] == an[0]
    ids = v[2:]; cls = np.array(an[2:])
    la = np.array([alat[x] for x in ids]); tc = np.array([tcls[x] for x in ids])
    is_high = cls == "high"
    nh = ~is_high
    s_inf = np.ones(len(ids))
    s_real = np.ones(len(ids))
    if nh.any():
        s_inf[nh] = pct(la[nh])
        s_real[nh] = pct(la[nh] + rng.normal(0, a.noise, nh.sum()))
    s_uninf = s_inf.copy()
    if nh.sum() > 1:
        s_uninf[nh] = rng.permutation(s_inf[nh])
    scores = {"informative": s_inf, "realistic": s_real, "uninformative": s_uninf}
    for k, s in scores.items():
        out[k] += [lines[i], lines[i + 1],
                   " ".join([v[0], "score:SIM"] + [f"{x:.4f}" for x in s])]
        stat[k].append(np.c_[s, tc, is_high])
    # synonymous-style: the low class only, scores drawn from the gene's medium
    # scores (informative arm) with replacement; a gene without both is skipped
    lo = cls == "low"; md = cls == "medium"
    if lo.any() and md.any():
        s_syn = rng.choice(s_inf[md], lo.sum(), replace=True)
        out["synonymous"] += [
            " ".join([v[0], "var"] + [x for x, m in zip(ids, lo) if m]),
            " ".join([v[0], "anno"] + ["low"] * int(lo.sum())),
            " ".join([v[0], "score:SIM"] + [f"{x:.4f}" for x in s_syn])]
        stat["synonymous"].append(np.c_[s_syn, tc[lo], np.zeros(lo.sum())])

for k in arms:
    with open(os.path.join(a.out, f"groupFile_{k}.txt"), "w") as fh:
        fh.write("\n".join(out[k]) + "\n")
    m = np.vstack(stat[k])
    nh = m[:, 2] == 0
    c_all = np.corrcoef(m[:, 0], m[:, 1])[0, 1] if len(m) > 2 else float("nan")
    c_nh = np.corrcoef(m[nh, 0], m[nh, 1])[0, 1] if nh.sum() > 2 else float("nan")
    print(f"{k:14s} genes={len(stat[k]):4d} variants={len(m):6d} "
          f"cor(score, true_class): all={c_all:.3f} medium/low={c_nh:.3f}")
