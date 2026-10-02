#!/usr/bin/env python3
"""Exact tail reference for a SAIGE spaBelow_/spaFail_ dump.

The dumpcheck's Monte Carlo already uses the right MODEL -- the carriers (nz)
as exact Bernoullis plus a Gaussian remainder over everyone else, which is the
fast kernel's own decomposition. Its problem is power: 2e6 draws bound the tail
at 5e-7 and the saddlepoint is claiming 1e-17, so it adjudicates nothing.

This replaces the sampling with an exact convolution in the SAME model:

    S = sum_{i in nz} gtilde_i x_i  +  Normal(rm, rs^2),   x_i ~ Bern(mu_i)

The discrete part is a convolution, not an enumeration -- O(n x grid), not 2^n.
Weights are rounded both DOWN and UP onto the grid, giving S_lo <= S <= S_hi
pointwise and therefore a RIGOROUS BRACKET on the tail rather than a point
estimate with an unquantified discretisation error.

Prints p-values only; no per-sample rows leave the function.
Usage:  python3 spa_exact_ref.py <dump.tsv> [...]
"""
import sys, math
import numpy as np


def pmf_grid(w, mu, step, mode):
    """pmf of sum(w_i x_i) on a lattice of `step`, weights rounded `mode`.
    Returns (offset_index, pmf) with value(j) = (offset + j) * step."""
    q = np.floor(w / step) if mode == "down" else np.ceil(w / step)
    q = q.astype(np.int64)
    lo = int(q[q < 0].sum())
    hi = int(q[q > 0].sum())
    size = hi - lo + 1
    pmf = np.zeros(size)
    pmf[-lo] = 1.0                      # S = 0 before any term
    for qi, m in zip(q, mu):
        if qi == 0:
            continue                    # contributes 0 either way
        nxt = pmf * (1.0 - m)
        if qi > 0:
            nxt[qi:] += pmf[:size - qi] * m
        else:
            nxt[:size + qi] += pmf[-qi:] * m
        pmf = nxt
    return lo, pmf


def tail(lo, pmf, step, q, rm, rs, side):
    """P(S_discrete + N(rm, rs^2) >= q)  (side=+1)  or  <= q  (side=-1)."""
    vals = (np.arange(len(pmf)) + lo) * step
    if rs <= 0:                # no Gaussian part: the remainder is a point mass
        ind = (vals + rm >= q) if side > 0 else (vals + rm <= q)
        return float(np.dot(pmf, ind.astype(float)))
    z = (q - vals - rm) / rs
    # survival of the normal part, accurate deep in the tail via erfc
    comp = 0.5 * np.array([math.erfc(zz / math.sqrt(2)) for zz in z]) if side > 0 \
      else 0.5 * np.array([math.erfc(-zz / math.sqrt(2)) for zz in z])
    return float(np.dot(pmf, comp))


def reference(gtilde, mu, nz, q, qinv, step=None, verbose=True):
    nz = nz.astype(bool)
    gNB, muNB = gtilde[nz], mu[nz]
    gNA, muNA = gtilde[~nz], mu[~nz]
    rm = float(np.sum(gNA * muNA))
    rs = float(np.sqrt(np.sum(muNA * (1 - muNA) * gNA ** 2)))
    if step is None:                      # ~1e5 lattice points across the range
        span = float(np.abs(gNB).sum())
        step = max(span / 2.0e5, 1e-12)
    out = {}
    for side, qq, name in ((+1, q, "P(S>=q)"), (-1, qinv, "P(S<=qinv)")):
        vals = []
        for mode in ("down", "up"):
            lo, pmf = pmf_grid(gNB, muNB, step, mode)
            vals.append(tail(lo, pmf, step, qq, rm, rs, side))
        out[name] = (min(vals), max(vals))
    if verbose:
        print(f"   exact-convolution reference (lattice {step:.3g}, "
              f"{int(nz.sum())} Bernoullis, remainder N({rm:.6g}, {rs:.4g}^2)):")
        for k, (a, b) in out.items():
            print(f"      {k:12s} in [{a:.4g}, {b:.4g}]")
        tot = (out["P(S>=q)"][0] + out["P(S<=qinv)"][0],
               out["P(S>=q)"][1] + out["P(S<=qinv)"][1])
        print(f"      {'two-sided':12s} in [{tot[0]:.4g}, {tot[1]:.4g}]")
    return out


def read_dump(path):
    hdr = {}
    with open(path) as f:
        l1 = f.readline().lstrip("# ").rstrip("\n").split("\t")
        hdr = dict(zip(l1[0::2], l1[1::2]))
        f.readline()                                   # line 2: kernel info
        cols = f.readline().rstrip("\n").split("\t")
    d = np.loadtxt(path, skiprows=3)
    ix = {c: i for i, c in enumerate(cols)}
    return (d[:, ix["gtilde"]], d[:, ix["mu"]], d[:, ix["nz"]],
            float(hdr["q"]), float(hdr["qinv"]))


if __name__ == "__main__":
    for path in sys.argv[1:]:
        g, mu, nz, q, qinv = read_dump(path)
        print(f"== {path.rsplit('/', 1)[-1]} | N={len(g)} nz={int(nz.sum())} "
              f"q={q:.6g} qinv={qinv:.6g}")
        reference(g, mu, nz, q, qinv)


def provenance(path):
    """Where a 15-order discrepancy between the counts and the p could live.

    Everything in reference() is DOWNSTREAM of gtilde: the saddlepoint, the
    convolution and the BETA/SE reconciliation all take the adjusted vector as
    given and agree because they ask the same question. This asks the question
    none of them do -- is gtilde consistent with the raw weighted burden g?

    Prints shares and ratios, no per-sample rows and no counts of people.
    """
    import numpy as np
    hdr = {}
    with open(path) as f:
        l1 = f.readline().lstrip("# ").rstrip("\n").split("\t")
        hdr = dict(zip(l1[0::2], l1[1::2]))
        f.readline()
        cols = f.readline().rstrip("\n").split("\t")
    d = np.loadtxt(path, skiprows=3)
    ix = {c: i for i, c in enumerate(cols)}
    gt, g, mu, y, nz = (d[:, ix["gtilde"]], d[:, ix["g"]], d[:, ix["mu"]],
                        d[:, ix["y"]], d[:, ix["nz"]].astype(bool))
    q = float(hdr["q"])

    def conc(v):
        a = np.sort(np.abs(v))[::-1]; s = a.sum()
        return (a[0] / s, a[:5].sum() / s, (s * s) / np.sum(a * a))

    print(f"-- {path.rsplit('/', 1)[-1]}")
    for name, v in (("gtilde", gt), ("g (raw)", g)):
        t1, t5, ne = conc(v[nz])
        print(f"   {name:9s} over nz: top1 {t1:6.1%}  top5 {t5:6.1%}  n_eff {ne:7.3f}  "
              f"max {np.abs(v[nz]).max():.4g}  sum|v| {np.abs(v[nz]).sum():.4g}")
    # the observation itself, against its own null mean -- model-free framing
    for name, v in (("gtilde", gt), ("g (raw)", g)):
        obs, exp = float(v @ y), float(v @ mu)
        sd = float(np.sqrt(np.sum(v * v * mu * (1 - mu))))
        print(f"   {name:9s} sum(v*y) {obs:12.6g}  E[sum(v*mu)] {exp:12.6g}  "
              f"sd {sd:10.4g}  z {(obs - exp) / sd:8.3f}")
    print(f"   q (header) {q:.6g}   gtilde: sum(v*(y-mu)) = {float(gt @ (y - mu)):.6g}")
