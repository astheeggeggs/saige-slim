#!/usr/bin/env python3
"""v6 against v4, one ancestry, from the gene-level summary tables (controls first).

Runs INSIDE a dsub job (the pyplot image) from the summary tool's outputs for the two
run tags, so nothing is read into or displayed from the notebook: PNGs and two small
gene-level tables are delivered. Everything here is a count of GENES or CELLS, never of
people, alleles or carriers.

Controls, stated before the numbers and printed first:
  1. synonymous: pooled QQ old vs new, and the fraction below 0.05 per trait -- both must
     sit on the null; inflation is the image or the model, not annotation.
  2. pLoF gene for gene: the pLoF rule did not change between v4 and v6, so the two
     runs' pLoF cells are the CROSS-IMAGE control. Spearman over shared (trait, gene),
     the count of cells moved by more than MOVE decades, and the hit counts.
     Expected movers: cells past 1e-5 (the saddlepoint tail and the exact-by-weight
     routing changed between the images, LEDGER #163-#176) and nothing else.
Contrasts (the result), each a scatter of -log10 p over shared (trait, gene) plus a
pooled QQ and the hit counts:
  3. splice: new pLoF;damaging_missense vs new pLoF;damaging_missense;spliceai2_gt040
  4. reannotation: old pLoF;damaging_missense vs new pLoF;damaging_missense;spliceai2_gt040
     (like for like: both carry a splice component, v4's inside damaging_missense at
     SpliceAI v1 >= 0.2) and vs new pLoF;damaging_missense (missense-only)
  5. AlphaMissense: new am_moderate vs new damaging_missense
  6. FlexRV: new p_flexrv (the Cauchy row) vs the new splice mask's SKAT-O p at CMP_MAF
Per trait, the summary table: hits past SIG per mask old and new, distinct genes, and the
synonymous fraction below 0.05. A p printed as 0 or below 1e-300 is clipped to 1e-300 for
the plots; the floor table (IN_FLOOR_NEW) lists those cells by their log10 column.

Env (dsub localizes IN_*, delocalizes OUT_DIR with --output-recursive):
  ANC            ancestry label
  IN_SUM_NEW     <anc>.flexrv_vs_skato.v6.tsv.gz     IN_CELLS_NEW  <anc>.skato_cells.v6.tsv.gz
  IN_SUM_OLD     <anc>.flexrv_vs_skato.v4.tsv.gz     IN_CELLS_OLD  <anc>.skato_cells.v4.tsv.gz
  IN_FLOOR_NEW   <anc>.floor.v6.tsv.gz (optional)
  OUT_DIR        directory for <anc>.controls.png, <anc>.contrasts.png, <anc>.v6_vs_v4.tsv,
                 <anc>.v6_vs_v4.txt
  CMP_MAF        the max_MAF cell compared (default 0.001); SIG (default 2.5e-6);
                 MOVE (default 1, decades); DPI (default 110)
Self-test: aou_v6_vs_v4_selftest.sh (fabricated tables, no AoU data).
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, beta

OLD = {"plof": "pLoF", "syn": "synonymous", "dmg": "pLoF;damaging_missense"}
NEW = {"plof": "pLoF", "syn": "synonymous", "dmg": "pLoF;damaging_missense",
       "spl": "pLoF;damaging_missense;spliceai2_gt040", "dm": "damaging_missense", "am": "am_moderate"}
FLOOR = 1e-300


def load(sum_path, cells_path, anc):
    S = pd.read_csv(sum_path, sep="\t", dtype=str); S = S[S.anc == anc]
    C = pd.read_csv(cells_path, sep="\t", dtype=str); C = C[C.anc == anc]
    return S, C


def col(C, kind, mask, maf):
    """The p-value series of one cell over (trait, gene), as floats; NaN where absent."""
    name = f"{kind}:{mask}@{maf}"
    if name not in C.columns:
        return pd.Series(np.nan, index=pd.MultiIndex.from_frame(C[["trait", "gene"]]), name=name)
    v = pd.to_numeric(C[name], errors="coerce")
    return pd.Series(v.values, index=pd.MultiIndex.from_frame(C[["trait", "gene"]]), name=name)


def neglog(p):
    return -np.log10(np.clip(p, FLOOR, 1.0))


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


def scatter_panel(ax, x, y, xl, yl, title, sig):
    j = x.index.intersection(y.index); xx, yy = x.loc[j], y.loc[j]
    ok = xx.notna() & yy.notna(); xx, yy = xx[ok], yy[ok]
    lx, ly = neglog(xx), neglog(yy)
    ax.scatter(lx, ly, s=4, color="#1f77b4", alpha=0.6)
    m = max(lx.max() if len(lx) else 1, ly.max() if len(ly) else 1)
    ax.plot([0, m], [0, m], color="0.3", lw=0.8)
    t = -np.log10(sig); ax.axhline(t, color="r", lw=0.5, ls=":"); ax.axvline(t, color="r", lw=0.5, ls=":")
    ax.set_xlabel(f"{xl}  $-\\log_{{10}} p$"); ax.set_ylabel(f"{yl}  $-\\log_{{10}} p$")
    ax.set_title(f"{title}\n{len(xx):,} shared (trait, gene); hits {int((xx < sig).sum())} vs {int((yy < sig).sum())}", fontsize=9)
    ax.grid(alpha=0.3)
    return xx, yy


def stats(x, y, sig, move):
    j = x.index.intersection(y.index); xx, yy = x.loc[j], y.loc[j]
    ok = xx.notna() & yy.notna(); xx, yy = xx[ok], yy[ok]
    if len(xx) < 3:
        return dict(n=len(xx), rho=np.nan, moved=0, hx=0, hy=0, both=0, only_x=0, only_y=0, genes_x=0, genes_y=0)
    lx, ly = neglog(xx), neglog(yy)
    hx, hy = xx < sig, yy < sig
    return dict(n=len(xx), rho=spearmanr(lx, ly).correlation, moved=int((np.abs(lx - ly) > move).sum()),
                hx=int(hx.sum()), hy=int(hy.sum()), both=int((hx & hy).sum()), only_x=int((hx & ~hy).sum()), only_y=int((~hx & hy).sum()),
                genes_x=xx[hx].index.get_level_values("gene").nunique(), genes_y=yy[hy].index.get_level_values("gene").nunique())


def main():
    anc = os.environ["ANC"]; out = os.environ["OUT_DIR"]; os.makedirs(out, exist_ok=True)
    maf = os.environ.get("CMP_MAF", "0.001"); sig = float(os.environ.get("SIG", "2.5e-6"))
    move = float(os.environ.get("MOVE", "1")); dpi = int(os.environ.get("DPI", "110"))
    Sn, Cn = load(os.environ["IN_SUM_NEW"], os.environ["IN_CELLS_NEW"], anc)
    So, Co = load(os.environ["IN_SUM_OLD"], os.environ["IN_CELLS_OLD"], anc)
    traits = sorted(set(Cn.trait) & set(Co.trait))
    lines = [f"{anc}: {len(traits)} traits in both runs ({Cn.trait.nunique()} new, {Co.trait.nunique()} old); cells compared at max_MAF {maf}; "
             f"hits past {sig:g}; a move is more than {move:g} decade(s)"]
    missing = [m for m in NEW.values() if f"burden:{m}@{maf}" not in Cn.columns]
    if missing:
        lines.append(f"WARNING: new cells table lacks masks {missing} at {maf}")
    old = {k: col(Co, "skato", m, maf) for k, m in OLD.items()}
    new = {k: col(Cn, "skato", m, maf) for k, m in NEW.items()}
    flex = pd.to_numeric(Sn.p_flexrv, errors="coerce"); flex.index = pd.MultiIndex.from_frame(Sn[["trait", "gene"]])

    # ---- controls -------------------------------------------------------------------
    lines.append("\n== CONTROLS (read these before the contrasts)")
    sy_o, sy_n = old["syn"].dropna(), new["syn"].dropna()
    lines.append(f"1. synonymous below 0.05: old {np.mean(sy_o < 0.05):.3f} (n={len(sy_o):,}), new {np.mean(sy_n < 0.05):.3f} (n={len(sy_n):,}); "
                 f"expected 0.050. Per-trait values in the TSV (syn_frac05_old/new).")
    st = stats(old["plof"], new["plof"], sig, move)
    lines.append(f"2. pLoF old vs new (cross-image control): n={st['n']:,} shared, Spearman {st['rho']:.4f}, moved >{move:g} decade: {st['moved']} "
                 f"(expected: only cells past 1e-5), hits old {st['hx']} new {st['hy']} both {st['both']} only-old {st['only_x']} only-new {st['only_y']}; "
                 f"distinct genes old {st['genes_x']} new {st['genes_y']}")
    # ---- contrasts ------------------------------------------------------------------
    lines.append("\n== CONTRASTS")
    pairs = [("3. splice label: new pLoF;dm vs new pLoF;dm;spliceai2_gt040", new["dmg"], new["spl"]),
             ("4a. reannotation, like for like: old pLoF;dm (v1 splice inside) vs new pLoF;dm;spliceai2_gt040", old["dmg"], new["spl"]),
             ("4b. reannotation, missense only: old pLoF;dm vs new pLoF;dm", old["dmg"], new["dmg"]),
             ("5. AlphaMissense: new damaging_missense vs new am_moderate", new["dm"], new["am"]),
             (f"6. FlexRV pooled vs new splice mask SKAT-O p at {maf}", new["spl"], flex)]
    for title, x, y in pairs:
        st = stats(x, y, sig, move)
        lines.append(f"{title}: n={st['n']:,}, Spearman {st['rho']:.3f}, moved {st['moved']}, hits {st['hx']} vs {st['hy']} "
                     f"(both {st['both']}, only-first {st['only_x']}, only-second {st['only_y']}); distinct genes {st['genes_x']} vs {st['genes_y']}")
    # ---- per-trait table ----------------------------------------------------------------
    rows = []
    for t in traits:
        r = dict(anc=anc, trait=t)
        for k, s in (("old_plof", old["plof"]), ("old_dmg", old["dmg"]), ("new_plof", new["plof"]), ("new_dmg", new["dmg"]),
                     ("new_spl", new["spl"]), ("new_dm", new["dm"]), ("new_am", new["am"])):
            v = s.loc[t] if t in s.index.get_level_values("trait") else pd.Series(dtype=float)
            r[f"hits_{k}"] = int((v.dropna() < sig).sum())
        v = flex.loc[t] if t in flex.index.get_level_values("trait") else pd.Series(dtype=float)
        r["hits_new_flexrv"] = int((v.dropna() < sig).sum())
        for k, s in (("syn_frac05_old", old["syn"]), ("syn_frac05_new", new["syn"])):
            v = s.loc[t].dropna() if t in s.index.get_level_values("trait") else pd.Series(dtype=float)
            r[k] = round(float(np.mean(v < 0.05)), 4) if len(v) else np.nan
        rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(os.path.join(out, f"{anc}.v6_vs_v4.tsv"), sep="\t", index=False)
    lines.append("\n== per trait (hits past threshold at the compared MAF; synonymous fraction below 0.05)")
    lines.append(T.to_string(index=False))
    if os.environ.get("IN_FLOOR_NEW") and os.path.exists(os.environ["IN_FLOOR_NEW"]):
        F = pd.read_csv(os.environ["IN_FLOOR_NEW"], sep="\t", dtype=str)
        lines.append(f"\n== floor: {len(F)} cell(s) past log10 p = -300 in the new run (by their log10 columns)")
        if len(F):
            lines.append(F.head(30).to_string(index=False))
    # ---- figures ------------------------------------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.5))
    qq_panel(ax[0], [("synonymous old (v4)", sy_o, "#7f7f7f"), ("synonymous new (v6)", sy_n, "#d62728")], f"{anc}: synonymous, pooled over traits, max_MAF {maf}")
    scatter_panel(ax[1], old["plof"], new["plof"], "pLoF old (v4)", "pLoF new (v6)", "pLoF gene for gene (cross-image control)", sig)
    qq_panel(ax[2], [("pLoF old (v4)", old["plof"].dropna(), "#7f7f7f"), ("pLoF new (v6)", new["plof"].dropna(), "#2ca02c")], "pLoF, pooled over traits")
    fig.suptitle(f"{anc.upper()}  controls: v6 against v4", y=0.995); fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(out, f"{anc}.controls.png"), dpi=dpi); plt.close(fig)
    fig, ax = plt.subplots(2, 3, figsize=(17, 10.5))
    scatter_panel(ax[0, 0], new["dmg"], new["spl"], "new pLoF;dm", "new pLoF;dm;spliceai2_gt040", "3. the splice label", sig)
    scatter_panel(ax[0, 1], old["dmg"], new["spl"], "old pLoF;dm", "new pLoF;dm;spliceai2_gt040", "4a. reannotation, like for like", sig)
    scatter_panel(ax[0, 2], old["dmg"], new["dmg"], "old pLoF;dm", "new pLoF;dm", "4b. reannotation, missense only", sig)
    scatter_panel(ax[1, 0], new["dm"], new["am"], "new damaging_missense", "new am_moderate", "5. AlphaMissense >= 0.906", sig)
    scatter_panel(ax[1, 1], new["spl"], flex, "new splice mask SKAT-O", "FlexRV pooled", "6. FlexRV against the splice mask", sig)
    qq_panel(ax[1, 2], [("old pLoF;dm", old["dmg"].dropna(), "#7f7f7f"), ("new pLoF;dm", new["dmg"].dropna(), "#1f77b4"),
                        ("new pLoF;dm;spliceai2_gt040", new["spl"].dropna(), "#d62728"), ("new am_moderate", new["am"].dropna(), "#ff7f0e"),
                        ("FlexRV pooled", flex.dropna(), "#9467bd")], "damaging masks, pooled over traits")
    fig.suptitle(f"{anc.upper()}  contrasts: v6 against v4 at max_MAF {maf}", y=0.995); fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(out, f"{anc}.contrasts.png"), dpi=dpi); plt.close(fig)
    txt = "\n".join(lines)
    with open(os.path.join(out, f"{anc}.v6_vs_v4.txt"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)
    print("done")


if __name__ == "__main__":
    main()
