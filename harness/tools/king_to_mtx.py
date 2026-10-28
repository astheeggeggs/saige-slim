#!/usr/bin/env python3
"""A sparse GRM for SAIGE from a KING-robust kinship table, and its comparison with an
existing SAIGE .mtx (LEDGER #183, 2026-10-28).

Why: amr's sparse GRM from createSparseGRM (5,000 random markers, standardised genotypes)
carries the admixture background in its off-diagonals -- 46M pairs above 0.05 -- which forced
the 0.125 cutoff, and a GRM thinned that hard is not positive definite (aou_sigma_pd_check.R).
KING-robust estimates kinship from within-pair genotype counts, so ancestry does not leak into
the off-diagonals; SAIGE accepts any sparse kinship in Matrix Market form, so this is an INPUT
change: no code path moves.

  convert   plink2 --make-king-table output (.kin0: #FID1 ID1 FID2 ID2 NSNP HETHET IBS0 KINSHIP,
            or the IID1/IID2 spelling) -> <out>.mtx on SAIGE's scale (off-diagonal = 2 x kinship,
            diagonal = 1, symmetric coordinate, lower triangle, 1-based) + <out>.mtx.sampleIDs.txt
            in the order of --samples (one IID per line, e.g. the .psam/.fam IIDs), keeping pairs
            at or above --cutoff on the GRM scale. An ID in the table but not in --samples is an
            error (the table was computed on these samples or it was not).
  compare   the KING-derived .mtx against an existing SAIGE .mtx on the SAME samples: (1) close
            relatives agree (pairs the existing GRM puts at or above 0.4: the new value's median
            ratio and spread -- the positive control; read it first); (2) pairs at or above each
            cutoff under each estimator (the admixture background is the difference); (3) the
            0.05..0.125 band: pairs both keep (real third-degree) against pairs only the old one
            keeps (background); (4) the diagonal of the old one. Counts of pairs below --mask are
            written "<mask" (real data); no sample ID is ever printed.
  selftest  fabricated data: three families planted, a background of small GRM-only entries,
            every claim checked.

Usage:
  king_to_mtx.py convert --king king.kin0 --samples ids.txt --out amr_king [--cutoff 0.05]
  king_to_mtx.py compare --new amr_king.mtx --old old.sparseGRM.mtx --report cmp.tsv [--mask 20]
  king_to_mtx.py selftest
"""
import argparse, gzip, io, os, sys, tempfile
import numpy as np


def openr(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_ids(p):
    with openr(p) as fh:
        ids = [l.strip().split()[0] for l in fh if l.strip()]
    if len(set(ids)) != len(ids):
        raise SystemExit("--samples has a repeated ID (SAIGE keys on IID alone)")
    return ids


def read_king(p):
    """(id1, id2, kinship) per row of a plink2 --make-king-table file."""
    with openr(p) as fh:
        head = fh.readline().lstrip("#").split()
        cols = {c.upper(): i for i, c in enumerate(head)}
        i1 = cols.get("ID1", cols.get("IID1")); i2 = cols.get("ID2", cols.get("IID2")); ik = cols.get("KINSHIP")
        if None in (i1, i2, ik):
            raise SystemExit(f"KING table header not recognised: {head}")
        for l in fh:
            f = l.split()
            if not f:
                continue
            yield f[i1], f[i2], float(f[ik])


def read_mtx(p):
    """A SAIGE .mtx: dict (i, j) with i >= j (0-based) -> value, plus n."""
    ent = {}
    with openr(p) as fh:
        head = fh.readline()
        if "symmetric" not in head:
            raise SystemExit("expected a symmetric coordinate .mtx")
        n = int(fh.readline().split()[0])
        for l in fh:
            f = l.split()
            if not f:
                continue
            i, j, x = int(f[0]) - 1, int(f[1]) - 1, float(f[2])
            if i < j:
                i, j = j, i
            ent[(i, j)] = x
    return n, ent


def write_mtx(path, n, ent):
    with open(path, "w") as fh:
        fh.write("%%MatrixMarket matrix coordinate real symmetric\n")
        fh.write(f"{n} {n} {len(ent)}\n")
        for (i, j) in sorted(ent):
            fh.write(f"{i + 1} {j + 1} {ent[(i, j)]!r}\n")


def convert(a):
    ids = read_ids(a.samples); pos = {s: k for k, s in enumerate(ids)}; n = len(ids)
    ent = {(k, k): 1.0 for k in range(n)}
    kept = dropped = unknown = 0
    for id1, id2, kin in read_king(a.king):
        if id1 not in pos or id2 not in pos:
            unknown += 1
            continue
        v = 2.0 * kin
        if v < a.cutoff:
            dropped += 1
            continue
        i, j = pos[id1], pos[id2]
        if i == j:
            continue
        if i < j:
            i, j = j, i
        ent[(i, j)] = v
        kept += 1
    if unknown:
        raise SystemExit(f"{unknown} KING pair(s) name a sample absent from --samples")
    write_mtx(a.out + ".mtx", n, ent)
    with open(a.out + ".mtx.sampleIDs.txt", "w") as fh:
        fh.write("\n".join(ids) + "\n")
    print(f"samples {n}; pairs kept at >= {a.cutoff} (GRM scale): {kept}; below: {dropped}")
    print(f"wrote {a.out}.mtx and {a.out}.mtx.sampleIDs.txt")


def mask(k, m):
    return f"<{m}" if 0 < k < m else str(k)


def compare(a):
    nN, new = read_mtx(a.new); nO, old = read_mtx(a.old)
    idsN = read_ids(a.new + ".sampleIDs.txt"); idsO = read_ids(a.old + ".sampleIDs.txt")
    if len(idsN) != nN or len(idsO) != nO:
        raise SystemExit("an .mtx and its sampleIDs file disagree on the sample count")
    # align on IID: a pair is compared on the samples both matrices carry
    posO = {s: k for k, s in enumerate(idsO)}
    common = [s for s in idsN if s in posO]
    rows = []
    def note(k, v):
        rows.append((k, str(v))); print(f"{k:48s} {v}")
    note("samples_new", nN); note("samples_old", nO); note("samples_common", len(common))
    mapN = {k: s for k, s in enumerate(idsN)}
    offN = {}
    for (i, j), x in new.items():
        if i != j:
            si, sj = mapN[i], mapN[j]
            if si in posO and sj in posO:
                a_, b_ = sorted((posO[si], posO[sj]), reverse=True)
                offN[(a_, b_)] = x
    offO = {k: x for k, x in old.items() if k[0] != k[1]}
    # (1) close relatives agree: the positive control
    close = [(offN.get(k, 0.0), x) for k, x in offO.items() if x >= 0.4]
    if close:
        r = np.array([nv / ov for nv, ov in close if nv > 0]); miss = sum(1 for nv, _ in close if nv == 0)
        note("old_pairs_at_or_above_0.4", mask(len(close), a.mask))
        note("  absent_from_new", mask(miss, a.mask))
        if len(r):
            note("  new_over_old_median", f"{np.median(r):.3f}")
            note("  new_over_old_5th_95th", f"{np.percentile(r, 5):.3f} .. {np.percentile(r, 95):.3f}")
    else:
        note("old_pairs_at_or_above_0.4", "0 (no first-degree pair in the old matrix on the common samples)")
    # (2) how many pairs each keeps at each cutoff
    for c in (0.05, 0.125, 0.25, 0.4):
        note(f"pairs_at_or_above_{c}_old", mask(sum(1 for x in offO.values() if x >= c), a.mask))
        note(f"pairs_at_or_above_{c}_new", mask(sum(1 for x in offN.values() if x >= c), a.mask))
    # (3) the 0.05..0.125 band: real third-degree against background
    band_old = {k for k, x in offO.items() if 0.05 <= x < 0.125}
    both = sum(1 for k in band_old if offN.get(k, 0.0) >= 0.05)
    note("band_0.05_0.125_old_pairs", mask(len(band_old), a.mask))
    note("  of which_new_also_at_or_above_0.05", mask(both, a.mask))
    note("  of which_new_below_0.05_(background)", mask(len(band_old) - both, a.mask))
    # (4) the old diagonal
    dO = np.array([x for (i, j), x in old.items() if i == j])
    if len(dO):
        note("old_diag_mean_min_max", f"{dO.mean():.4f} {dO.min():.4f} {dO.max():.4f}")
    note("old_diag_entries_missing", mask(nO - len(dO), a.mask))
    with open(a.report, "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in rows:
            fh.write(f"{k.strip()}\t{v}\n")
    print("done")


def selftest(_):
    rng = np.random.default_rng(5); n = 300
    ids = [f"s{k}" for k in range(n)]
    T = tempfile.mkdtemp()
    with open(f"{T}/ids.txt", "w") as fh:
        fh.write("\n".join(ids) + "\n")
    # truth: three families (a trio-ish star at 0.5, a sib pair at 0.5, two cousins at 0.125)
    fam = {(1, 0): 0.5, (2, 0): 0.5, (11, 10): 0.5, (21, 20): 0.125}
    # the KING table: true pairs plus noise, kinship = GRM/2
    with open(f"{T}/king.kin0", "w") as fh:
        fh.write("#FID1\tID1\tFID2\tID2\tNSNP\tHETHET\tIBS0\tKINSHIP\n")
        for (i, j), v in fam.items():
            fh.write(f"0\t{ids[i]}\t0\t{ids[j]}\t100000\t0\t0\t{v / 2 + rng.normal(0, 0.002):.5f}\n")
        for _ in range(50):   # unrelated pairs KING reports near zero (some slightly negative)
            i, j = rng.choice(n, 2, replace=False)
            fh.write(f"0\t{ids[i]}\t0\t{ids[j]}\t100000\t0\t0\t{rng.normal(0, 0.005):.5f}\n")
    # the old GRM: the families plus an admixture background of 400 entries at 0.05..0.12
    old = {(k, k): 1.0 + rng.normal(0, 0.01) for k in range(n)}
    for (i, j), v in fam.items():
        old[(i, j)] = v + rng.normal(0, 0.01)
    bg = 0
    while bg < 400:
        i, j = sorted(rng.choice(n, 2, replace=False), reverse=True)
        if (i, j) not in old:
            old[(i, j)] = rng.uniform(0.05, 0.12); bg += 1
    write_mtx(f"{T}/old.mtx", n, old)
    with open(f"{T}/old.mtx.sampleIDs.txt", "w") as fh:
        fh.write("\n".join(ids) + "\n")
    convert(argparse.Namespace(king=f"{T}/king.kin0", samples=f"{T}/ids.txt", out=f"{T}/new", cutoff=0.05))
    nN, new = read_mtx(f"{T}/new.mtx")
    off = {k: x for k, x in new.items() if k[0] != k[1]}
    assert nN == n and len(off) == 4, f"expected the 4 planted pairs at >= 0.05, got {len(off)}"
    assert all(abs(off[k] - v) < 0.02 for k, v in fam.items()), "planted pairs not on the GRM scale"
    assert all(new[(k, k)] == 1.0 for k in range(n)), "diagonal must be 1"
    buf = io.StringIO(); so = sys.stdout; sys.stdout = buf
    try:
        compare(argparse.Namespace(new=f"{T}/new.mtx", old=f"{T}/old.mtx", report=f"{T}/cmp.tsv", mask=0))
    finally:
        sys.stdout = so
    rep = dict(l.rstrip("\n").split("\t") for l in open(f"{T}/cmp.tsv") if not l.startswith("key"))
    assert rep["old_pairs_at_or_above_0.4"] == "3", rep
    assert rep["absent_from_new"] == "0"
    assert 0.95 < float(rep["new_over_old_median"]) < 1.05, rep["new_over_old_median"]
    assert rep["pairs_at_or_above_0.05_old"] == str(400 + 4) and rep["pairs_at_or_above_0.05_new"] == "4"
    assert rep["pairs_at_or_above_0.25_old"] == "3" and rep["pairs_at_or_above_0.25_new"] == "3"
    assert rep["band_0.05_0.125_old_pairs"] == "400", rep["band_0.05_0.125_old_pairs"]
    assert rep["of which_new_also_at_or_above_0.05"] == "0" and rep["of which_new_below_0.05_(background)"] == "400"
    assert rep["old_diag_entries_missing"] == "0"
    # masking: a count below the mask is written "<mask"
    buf = io.StringIO(); sys.stdout = buf
    try:
        compare(argparse.Namespace(new=f"{T}/new.mtx", old=f"{T}/old.mtx", report=f"{T}/cmp2.tsv", mask=20))
    finally:
        sys.stdout = so
    rep2 = dict(l.rstrip("\n").split("\t") for l in open(f"{T}/cmp2.tsv") if not l.startswith("key"))
    assert rep2["old_pairs_at_or_above_0.4"] == "<20" and rep2["pairs_at_or_above_0.05_old"] == "404", rep2
    print("selftest OK: convert (4 planted pairs on the GRM scale, diagonal 1), compare (positive control "
          "ratio ~1, 400 background pairs only in the old matrix, the band split, masking)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("convert"); c.add_argument("--king", required=True); c.add_argument("--samples", required=True)
    c.add_argument("--out", required=True); c.add_argument("--cutoff", type=float, default=0.05); c.set_defaults(fn=convert)
    k = sub.add_parser("compare"); k.add_argument("--new", required=True); k.add_argument("--old", required=True)
    k.add_argument("--report", required=True); k.add_argument("--mask", type=int, default=20); k.set_defaults(fn=compare)
    s = sub.add_parser("selftest"); s.set_defaults(fn=selftest)
    a = ap.parse_args(); a.fn(a)


if __name__ == "__main__":
    main()
