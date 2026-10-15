#!/usr/bin/env python3
"""Per-(ancestry, trait) QQ grids from the gene-level summary tables, rendered
INSIDE a dsub job (image gcr.io/<project>/pyplot) so that nothing is read into
or displayed from the notebook: PNGs only are delivered.

Six panels, one per SKAT-O annotation mask, each the QQ of the burden p over
genes at the three max MAFs with the Beta order-statistic 95% band; a seventh
with the six masks at max MAF 0.001 plus the FlexRV pooled p. The bulk of each
series is thinned to ~1,500 evenly spaced ranks below expected 2; every point
above is kept, so the tails -- the figure -- are exact.

Env (dsub localizes IN_*, delocalizes OUT_DIR with --output-recursive):
  ANC          ancestry label for titles and file names
  IN_SUMMARY   <anc>.flexrv_vs_skato.<tag>.tsv.gz   (p_flexrv per trait x gene)
  IN_CELLS     <anc>.skato_cells.<tag>.tsv.gz       (burden:<mask>@<maf> columns)
  OUT_DIR      directory to write <anc>_<trait>.png into
  DPI          default 110
Self-test: aou_flexrv_qq_selftest.sh (fabricated tables, no AoU data).
"""
import os, sys, gzip
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import beta

# The masks are whatever the cells table holds (its columns are `burden:<mask>@<maf>`),
# drawn in this preferred order with anything unlisted appended alphabetically. v4's
# six masks first, then v6's (2026-10-15: the seven-mask set with SpliceAI2 at 0.40 and
# the am_moderate arm, which the summary tool folds into the same table).
MASK_PREF = ["pLoF", "damaging_missense", "other_missense", "synonymous",
             "pLoF;damaging_missense", "pLoF;damaging_missense;spliceai2_gt040",
             "pLoF;damaging_missense;other_missense;synonymous",
             "pLoF;damaging_missense;spliceai2_gt040;other_missense;synonymous", "am_moderate"]
def masks_in(C):
    have = {c[len("burden:"):].rsplit("@", 1)[0] for c in C.columns if c.startswith("burden:")}
    return [m for m in MASK_PREF if m in have] + sorted(have - set(MASK_PREF))
FLEX_MAF = 0.001
TRAITS_BINARY = {"AMD", "Asth", "AFib", "BenCervUterNeo", "BenIntNeo", "BreastCanc", "CervCanc",
                 "COPD", "CRF", "ColonRectCanc", "CAD", "EFRMB", "FemInf", "Gout", "HF", "HTN", "IBD",
                 "IFHern", "ILDSarc", "MatHem", "NonRheuValv", "Pancreat", "PeptUlcer", "PAD",
                 "Psori", "RheumHeaDis", "RheumArth", "Stroke", "T2Diab", "Urolith", "VaricVeins", "VTE"}

def nice(s): return s.replace("_", " ")

def qq_points(p, keep_top=2.0, bulk=1500):
    p = np.sort(np.asarray(p, float)); p = p[np.isfinite(p) & (p > 0)]
    n = len(p)
    if n == 0: return np.array([]), np.array([]), 0
    k = np.arange(1, n + 1); e = -np.log10((k - 0.5) / n); o = -np.log10(p)
    top = e >= keep_top
    idx = np.where(top)[0]; rest = np.where(~top)[0]
    if len(rest) > bulk: rest = rest[np.linspace(0, len(rest) - 1, bulk).astype(int)]
    sel = np.sort(np.concatenate([idx, rest]))
    return e[sel], o[sel], n

def null_band(n, bulk=400):
    k = np.unique(np.concatenate([np.arange(1, min(n, 200) + 1), np.linspace(1, n, bulk).astype(int)]))
    e = -np.log10((k - 0.5) / n)
    lo = -np.log10(beta.ppf(0.975, k, n + 1 - k)); hi = -np.log10(beta.ppf(0.025, k, n + 1 - k))
    return e, lo, hi

def qq_grid(anc, trait, c, s, masks):
    """c: this trait's rows of the cells table; s: this trait's rows of the summary; masks: the panel list."""
    cells = {}
    for col in [x for x in c.columns if x.startswith("burden:")]:
        g, m = col[len("burden:"):].rsplit("@", 1); cells.setdefault(g, {})[float(m)] = col
    mafs = sorted({m for d in cells.values() for m in d})
    blues = plt.cm.Blues(np.linspace(0.35, 0.95, max(len(mafs), 1)))
    nrow = max(2, -(-(len(masks) + 1) // 4))   # the mask panels plus the FlexRV panel, four across
    fig = plt.figure(figsize=(18, 4.25 * nrow))
    axes = fig.subplots(nrow, 4)       # independent axes: a p of 1e-78 (BRCA1) must not set every panel's scale, nor be clipped
    fig.suptitle(f"{anc.upper()}  {trait}   burden $p$-value per annotation mask and max MAF (SAIGE-GENE+ cells), and FlexRV pooled vs all masks at maxMAF {FLEX_MAF:g}", y=0.995)
    for ax, g in zip(axes.flat[:len(masks)], masks):
        ax.set_title(nice(g).replace(";", ";\n") if g.count(";") > 1 else nice(g), fontsize=10)
        nmax = 0
        for m, col in zip(mafs, blues):
            colname = cells.get(g, {}).get(m)
            if colname is None: continue
            e, o, n = qq_points(pd.to_numeric(c[colname], errors="coerce").dropna())
            if n == 0: continue
            ax.scatter(e, o, s=4, color=col, label=f"maxMAF {m:g}  (n={n:,})"); nmax = max(nmax, n)
        if nmax:
            e, lo, hi = null_band(nmax); ax.fill_between(e, lo, hi, color="grey", alpha=0.25, lw=0, label="95% null band")
            ax.plot([0, e.max()], [0, e.max()], color="0.3", lw=0.8)
        ax.grid(alpha=0.3); ax.legend(fontsize=7, loc="upper left")
    ax = axes.flat[len(masks)]; ax.set_title(f"FlexRV pooled vs burden per mask at maxMAF {FLEX_MAF:g}", fontsize=10)
    palette = ["#08306b", "#2ca02c", "#ff7f0e", "#7f7f7f", "#6a3d9a", "#8c564b", "#17becf", "#bcbd22", "#e377c2"]
    nmax = 0
    for g, col in zip(masks, palette):
        colname = cells.get(g, {}).get(FLEX_MAF)
        if colname is None: continue
        e, o, n = qq_points(pd.to_numeric(c[colname], errors="coerce").dropna())
        if n: ax.scatter(e, o, s=4, color=col, label=(nice(g) if len(g) < 30 else nice(g)[:28] + "...")); nmax = max(nmax, n)
    e, o, n = qq_points(pd.to_numeric(s.p_flexrv, errors="coerce").dropna())
    if n: ax.scatter(e, o, s=5, color="#d62728", label=f"FlexRV pooled (n={n:,})"); nmax = max(nmax, n)
    if nmax:
        e, lo, hi = null_band(nmax); ax.fill_between(e, lo, hi, color="grey", alpha=0.25, lw=0, label="95% null band")
        ax.plot([0, e.max()], [0, e.max()], color="0.3", lw=0.8)
    ax.grid(alpha=0.3); ax.legend(fontsize=7, loc="upper left")
    for ax in axes.flat[len(masks) + 1:]: ax.axis("off")
    for ax in axes.flat[:len(masks) + 1]:
        ax.set_xlim(left=0); ax.set_ylim(bottom=0)
        ax.set_xlabel("expected $-\\log_{10} p$"); ax.set_ylabel("observed $-\\log_{10} p$")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig

def main():
    anc = os.environ["ANC"]; out = os.environ["OUT_DIR"]; dpi = int(os.environ.get("DPI", "110"))
    os.makedirs(out, exist_ok=True)
    S = pd.read_csv(os.environ["IN_SUMMARY"], sep="\t", dtype=str)
    C = pd.read_csv(os.environ["IN_CELLS"], sep="\t", dtype=str)
    S = S[S.anc == anc]
    traits = sorted(set(C.trait) & set(S.trait), key=lambda t: (t not in TRAITS_BINARY, t))
    masks = masks_in(C)
    print(f"{anc}: {len(traits)} traits with both tables; cells columns: {sum(x.startswith('burden:') for x in C.columns)}")
    print(f"masks: {len(masks)} -> " + ", ".join(masks))
    for t in traits:
        fig = qq_grid(anc, t, C[C.trait == t], S[S.trait == t], masks)
        fig.savefig(os.path.join(out, f"{anc}_{t}.png"), dpi=dpi); plt.close(fig)
    print(f"wrote {len(traits)} figures to {out}")
    print("done")

if __name__ == "__main__":
    main()
