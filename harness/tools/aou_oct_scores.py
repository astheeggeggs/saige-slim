#!/usr/bin/env python3
"""The October release's FlexRV scores against each other and against their null, one
ancestry, from the gene-level summary tables (controls first).

Runs INSIDE a dsub job (the pyplot image) from aou_flexrv_summary.sh's outputs for the
three FlexRV arms of the October fleet, so nothing is read into or displayed from the
notebook: one PNG, one text report and one per-trait table are delivered. Everything here
is a count of GENES or CELLS, never of people, alleles or carriers.

Why: aou_v6_vs_v4.py takes ONE FlexRV table, so with the AMPCT arm given it never sees
SCALLION, and the per-trait QQ grids draw one FlexRV series each. The release's own
question -- do the two percentile scores find the same genes, and does either find more
than a random score would -- had no figure.

Controls, printed first:
  1. RAND on synonymous (the machinery's null, 19.9): the observed p at the expected
     1e-1 .. 1e-5 points and the hits past SIG. Must sit on the null; the body (fraction
     below 0.05) is the Cauchy's and is reported but not read. A RAND trait whose own tail
     carries a gene below RAND_LEAK is LEAKY -- the unweighted-score sets are a real
     synonymous burden, and a polygenic trait (Height, 2026-10-28, every ancestry) shows
     through them whatever the score -- so it is named, with its gene count below RAND_LEAK,
     and left out of the pooled null line and the null panels; the all-traits pool is
     printed beside it for reference. Read the null from the clean pool.
  2. the SKAT-O synonymous burden at CMP_MAF from the cells table (optional, IN_CELLS):
     the model's null beside the machinery's.
Contrast: AMPCT against SCALLION over shared (trait, gene): Spearman of -log10 p, cells
moved by more than MOVE decades, hits past SIG on each side, both, only-one; then each
score's pooled tail quantiles beside RAND's. RAND ran on a DIFFERENT group (synonymous),
so it is compared as a pooled tail only, never gene for gene.

Env (dsub localizes IN_*, delocalizes OUT_DIR with --output-recursive):
  ANC            ancestry label
  IN_AMPCT       <anc>.flexrv_vs_skato.oct_fxa.tsv.gz      (p_flexrv per trait x gene)
  IN_SCALLION    <anc>.flexrv_vs_skato.oct_fxs.tsv.gz
  IN_RAND        <anc>.flexrv_vs_skato.oct_rand.tsv.gz     (optional; fewer traits)
  IN_CELLS       <anc>.skato_cells.oct_fxa.tsv.gz          (optional; synonymous burden)
  OUT_DIR        directory for <anc>.scores.png, <anc>.scores.tsv, <anc>.scores.txt
  CMP_MAF        the synonymous cell read from IN_CELLS (default 0.001); SIG (default
                 2.5e-6); MOVE (default 1, decades); DPI (default 110); RAND_LEAK (default
                 1e-8): a RAND trait with any gene below it is leaky and leaves the pooled null
Self-test: aou_oct_scores_selftest.sh (fabricated tables, no AoU data).
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, beta

FLOOR = 1e-300
EXPECTED = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]


def load_flex(path, anc):
    S = pd.read_csv(path, sep="\t", dtype=str); S = S[S.anc == anc]
    v = pd.to_numeric(S.p_flexrv, errors="coerce")
    return pd.Series(v.values, index=pd.MultiIndex.from_frame(S[["trait", "gene"]]), name=os.path.basename(path))


def load_cell(path, anc, mask, maf):
    C = pd.read_csv(path, sep="\t", dtype=str); C = C[C.anc == anc]
    name = f"burden:{mask}@{maf}"
    if name not in C.columns:
        return None
    v = pd.to_numeric(C[name], errors="coerce")
    return pd.Series(v.values, index=pd.MultiIndex.from_frame(C[["trait", "gene"]]), name=name)


def neglog(p):
    return -np.log10(np.clip(p, FLOOR, 1.0))


def tail(p, sig):
    """Observed p at the expected 1e-k points (the k-th quantile of n p-values), the hits
    past sig and the fraction below 0.05; the quantiles at points the n cannot resolve are NaN."""
    p = np.sort(np.asarray(p, float)); p = p[np.isfinite(p)]; n = len(p)
    q = {}
    for e in EXPECTED:
        k = int(np.ceil(n * e))
        q[e] = p[k - 1] if k >= 1 and n else np.nan
    return dict(n=n, q=q, hits=int((p < sig).sum()) if n else 0, f05=float(np.mean(p < 0.05)) if n else np.nan)


def tail_line(label, t):
    qs = " / ".join(f"{t['q'][e]:.1e}" if np.isfinite(t['q'][e]) else "NA" for e in EXPECTED)
    return f"{label}: n={t['n']:,}; observed p at expected 1e-1/1e-2/1e-3/1e-4/1e-5: {qs}; hits {t['hits']}; below 0.05 {t['f05']:.3f}"


def qq_points(p, keep_top=2.0, bulk=1500):
    p = np.sort(np.asarray(p, float)); p = p[np.isfinite(p)]
    p = np.clip(p, FLOOR, 1.0); n = len(p)
    if n == 0:
        return np.array([]), np.array([]), 0
    k = np.arange(1, n + 1); e = -np.log10((k - 0.5) / n); o = -np.log10(p)
    top = e >= keep_top
    idx = np.where(top)[0]; rest = np.where(~top)[0]
    if len(rest) > bulk:
        rest = rest[np.linspace(0, len(rest) - 1, bulk).astype(int)]
    sel = np.sort(np.concatenate([idx, rest]))
    return e[sel], o[sel], n


def null_band(n, bulk=400):
    k = np.unique(np.concatenate([np.arange(1, min(n, 200) + 1), np.linspace(1, n, bulk).astype(int)]))
    e = -np.log10((k - 0.5) / n)
    return e, -np.log10(beta.ppf(0.975, k, n + 1 - k)), -np.log10(beta.ppf(0.025, k, n + 1 - k))


def qq_panel(ax, series, title):
    nmax = 0
    for label, p, c in series:
        e, o, n = qq_points(p)
        if n:
            ax.scatter(e, o, s=4, color=c, label=f"{label} (n={n:,})"); nmax = max(nmax, n)
    if nmax:
        e, lo, hi = null_band(nmax); ax.fill_between(e, lo, hi, color="grey", alpha=0.25, lw=0)
        ax.plot([0, e.max()], [0, e.max()], color="0.3", lw=0.8)
    ax.set_title(title, fontsize=10); ax.grid(alpha=0.3); ax.legend(fontsize=7, loc="upper left")
    ax.set_xlabel("expected $-\\log_{10} p$"); ax.set_ylabel("observed $-\\log_{10} p$")


def shared(x, y):
    j = x.index.intersection(y.index); xx, yy = x.loc[j], y.loc[j]
    ok = xx.notna() & yy.notna()
    return xx[ok], yy[ok]


def stats(x, y, sig, move):
    xx, yy = shared(x, y)
    if len(xx) < 3:
        return dict(n=len(xx), rho=np.nan, moved=0, hx=0, hy=0, both=0, only_x=0, only_y=0, genes_x=0, genes_y=0)
    lx, ly = neglog(xx), neglog(yy)
    hx, hy = xx < sig, yy < sig
    return dict(n=len(xx), rho=spearmanr(lx, ly).correlation, moved=int((np.abs(lx - ly) > move).sum()),
                hx=int(hx.sum()), hy=int(hy.sum()), both=int((hx & hy).sum()), only_x=int((hx & ~hy).sum()), only_y=int((~hx & hy).sum()),
                genes_x=xx[hx].index.get_level_values("gene").nunique(), genes_y=yy[hy].index.get_level_values("gene").nunique())


def scatter_panel(ax, x, y, xl, yl, title, sig):
    xx, yy = shared(x, y)
    lx, ly = neglog(xx), neglog(yy)
    ax.scatter(lx, ly, s=4, color="#1f77b4", alpha=0.6)
    m = max(lx.max() if len(lx) else 1, ly.max() if len(ly) else 1)
    ax.plot([0, m], [0, m], color="0.3", lw=0.8)
    t = -np.log10(sig); ax.axhline(t, color="r", lw=0.5, ls=":"); ax.axvline(t, color="r", lw=0.5, ls=":")
    ax.set_xlabel(f"{xl}  $-\\log_{{10}} p$"); ax.set_ylabel(f"{yl}  $-\\log_{{10}} p$")
    ax.set_title(f"{title}\n{len(xx):,} shared (trait, gene); hits {int((xx < sig).sum())} vs {int((yy < sig).sum())}", fontsize=9)
    ax.grid(alpha=0.3)


def main():
    anc = os.environ["ANC"]; out = os.environ["OUT_DIR"]; os.makedirs(out, exist_ok=True)
    maf = os.environ.get("CMP_MAF", "0.001"); sig = float(os.environ.get("SIG", "2.5e-6"))
    move = float(os.environ.get("MOVE", "1")); dpi = int(os.environ.get("DPI", "110"))
    leak = float(os.environ.get("RAND_LEAK", "1e-8"))
    am = load_flex(os.environ["IN_AMPCT"], anc); sc = load_flex(os.environ["IN_SCALLION"], anc)
    rand = load_flex(os.environ["IN_RAND"], anc) if os.environ.get("IN_RAND") and os.path.exists(os.environ["IN_RAND"]) else None
    syn = load_cell(os.environ["IN_CELLS"], anc, "synonymous", maf) if os.environ.get("IN_CELLS") and os.path.exists(os.environ["IN_CELLS"]) else None
    traits = sorted(set(am.index.get_level_values("trait")) & set(sc.index.get_level_values("trait")))
    # a RAND trait with a gene below RAND_LEAK carries signal through the unweighted sets: name it, pool without it
    leaky = {}
    if rand is not None:
        for t in sorted(set(rand.index.get_level_values("trait"))):
            v = rand.loc[t].dropna()
            if len(v) and v.min() < leak:
                leaky[t] = (float(v.min()), int((v < leak).sum()))
    rand_clean = rand[~rand.index.get_level_values("trait").isin(leaky)] if rand is not None else None
    rand_label = "RAND on synonymous" + (f" (excl. {len(leaky)} leaky)" if leaky else "")
    lines = [f"{anc}: {len(traits)} traits with both scores (AMPCT {am.index.get_level_values('trait').nunique()}, "
             f"SCALLION {sc.index.get_level_values('trait').nunique()}, RAND {rand.index.get_level_values('trait').nunique() if rand is not None else 0}); "
             f"hits past {sig:g}; a move is more than {move:g} decade(s)"]
    # ---- controls -------------------------------------------------------------------
    lines.append("\n== CONTROLS (read these before the contrast; the tail quantiles and the hits, not the fraction below 0.05)")
    if rand is not None:
        if leaky:
            lines.append(f"1. RAND traits with a gene below {leak:g} (leaky: signal through the unweighted sets; left out of the pooled null): "
                         + ", ".join(f"{t} (min p {m:.1e}, {k} gene{'s' if k != 1 else ''})" for t, (m, k) in leaky.items()))
            lines.append("   " + tail_line(f"RAND on synonymous, pooled over the {rand_clean.index.get_level_values('trait').nunique()} clean traits (the machinery's null: READ THIS ONE)", tail(rand_clean.dropna(), sig)))
            lines.append("   " + tail_line("RAND on synonymous, pooled over ALL its traits (for reference)", tail(rand.dropna(), sig)))
        else:
            lines.append(f"1. no RAND trait has a gene below {leak:g}")
            lines.append("   " + tail_line("RAND on synonymous, pooled over its traits (the machinery's null)", tail(rand.dropna(), sig)))
    else:
        lines.append("1. RAND table not given: the machinery's null is not shown")
    if syn is not None:
        lines.append("2. " + tail_line(f"SKAT-O synonymous burden at {maf}, pooled over traits (the model's null)", tail(syn.dropna(), sig)))
    else:
        lines.append("2. cells table not given or lacks the synonymous cell: the model's null is not shown")
    # ---- contrast -------------------------------------------------------------------
    lines.append("\n== CONTRAST: AMPCT against SCALLION (gene for gene) and each against RAND (pooled tail only: RAND ran on another group)")
    st = stats(am, sc, sig, move)
    lines.append(f"3. AMPCT vs SCALLION: n={st['n']:,} shared, Spearman {st['rho']:.3f}, moved >{move:g} decade: {st['moved']}, "
                 f"hits AMPCT {st['hx']} SCALLION {st['hy']} both {st['both']} only-AMPCT {st['only_x']} only-SCALLION {st['only_y']}; "
                 f"distinct genes AMPCT {st['genes_x']} SCALLION {st['genes_y']}")
    lines.append("4. " + tail_line("AMPCT pooled over traits", tail(am.dropna(), sig)))
    lines.append("5. " + tail_line("SCALLION pooled over traits", tail(sc.dropna(), sig)))
    # ---- per-trait table ----------------------------------------------------------------
    rows = []
    for t in traits:
        a = am.loc[t].dropna(); s = sc.loc[t].dropna()
        r = dict(anc=anc, trait=t, n_genes=len(a.index.intersection(s.index)))
        ta, ts = tail(a, sig), tail(s, sig)
        r.update(hits_ampct=ta["hits"], hits_scallion=ts["hits"], hits_both=int(len(a[a < sig].index.intersection(s[s < sig].index))),
                 q1e3_ampct=ta["q"][1e-3], q1e4_ampct=ta["q"][1e-4], q1e3_scallion=ts["q"][1e-3], q1e4_scallion=ts["q"][1e-4])
        if rand is not None and t in rand.index.get_level_values("trait"):
            v = rand.loc[t].dropna(); tr = tail(v, sig)
            r.update(hits_rand=tr["hits"], q1e3_rand=tr["q"][1e-3], q1e4_rand=tr["q"][1e-4],
                     rand_min_p=float(v.min()) if len(v) else np.nan, rand_leaky=t in leaky)
        else:
            r.update(hits_rand=np.nan, q1e3_rand=np.nan, q1e4_rand=np.nan, rand_min_p=np.nan, rand_leaky=np.nan)
        if syn is not None and t in syn.index.get_level_values("trait"):
            r["syn_frac05_skato"] = round(float(np.mean(syn.loc[t].dropna() < 0.05)), 4)
        rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(out, f"{anc}.scores.tsv"), sep="\t", index=False, float_format="%.3g")
    lines.append("\n== per trait (hits past threshold; observed p at the expected 1e-3 and 1e-4 points)")
    lines.append(T.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    # ---- figure ------------------------------------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.5))
    series = [(rand_label, rand_clean.dropna(), "#7f7f7f")] if rand is not None else []
    if syn is not None:
        series.append((f"SKAT-O synonymous @ {maf}", syn.dropna(), "#bcbd22"))
    qq_panel(ax[0], series, f"{anc}: the nulls, pooled over traits")
    scatter_panel(ax[1], am, sc, "AMPCT", "SCALLION", "AMPCT against SCALLION, p_FlexRV", sig)
    qq_panel(ax[2], ([(rand_label, rand_clean.dropna(), "#7f7f7f")] if rand is not None else []) +
             [("AMPCT", am.dropna(), "#1f77b4"), ("SCALLION", sc.dropna(), "#d62728")], "the two scores and the null, pooled over traits")
    fig.suptitle(f"{anc.upper()}  October FlexRV scores", y=0.995); fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(out, f"{anc}.scores.png"), dpi=dpi); plt.close(fig)
    txt = "\n".join(lines)
    with open(os.path.join(out, f"{anc}.scores.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)
    print("done")


if __name__ == "__main__":
    main()
