#!/usr/bin/env python3
"""FlexRV power fixture, second generation (2026-09-14): phenotypes only.

The 2026-09-10 fixture put 27 of its 62 causal genes' entire visible signal
below MAC 20 (handoff-2026-09-13-flexrv-sweep-rerun.md, 6c-6e), so the
oracle burden could detect only a handful of genes and no score arm could be
separated from another. This script leaves every genotype, GRM, PC, group
file and score line untouched and writes NEW phenotype columns whose causal
signal sits where a burden test can resolve it:

  * causal genes are drawn among genes whose VISIBLE cell (anno high or
    medium, MAF <= --maxmaf -- the cell every study run tests) carries at
    least --min-totmac minor alleles;
  * within a causal gene P(causal) is a steep logistic in the within-gene
    percentile of a_lat, the latent the informative score is the rank of
    (make_scores.py), with the `high` class at the top of it -- so the
    informative score IS informative about which variants carry effect, the
    realistic score (a_lat + noise) partly, and the permuted score not at all;
  * every causal variant of a gene has the SAME beta and the same sign
    (directional: the burden is the test under study), so the causal-only
    burden is the gene's oracle and the per-gene non-centrality is
    N * beta^2 * sum 2p(1-p) over the causal variants, drawn log-uniform on
    --ladder and delivered exactly on the realised burden;
  * the background is the fixture's own y_null (pedigree-transmitted
    polygenic + covariates + noise, variance 1), scaled by sqrt(1 - h2), so
    the new traits keep the GRM structure the null models were built on;
  * the binary arm is at --prevalence (default 0.10, not the fixture's 0.01:
    at K = 0.01 a burden test on 500 cases has 0.072 of the quantitative
    non-centrality and the equal-NCP liability h2 does not fit the budget)
    with its OWN liability, each gene scaled so that the EXACT
    liability-threshold non-centrality of its burden equals the quantitative
    one (the linearised constant saturates above ~0.6 sd per allele);
    y_null_bin<suffix> is y_null thresholded at the same K (the 2026-09-08
    recipe; permuting is ruled out, LEDGER #104).

Prediction first: with the ACTUAL score lines of the three arms read from
--scores, the script prints, before writing anything, the expected number of
Bonferroni detections (alpha 0.05 / n_genes, 1-df non-central chi-square) for
the unweighted burden, the causal-only oracle and each arm's score-weighted
burden at several exponents. The sweep is then compared against that table.

Writes, beside the fixture (nothing existing is modified):
  pheno_pcs<suffix>.txt   every column of pheno_pcs.txt plus
                          y_quant<suffix> y_binary<suffix> y_null_bin<suffix>
  truth_genes<suffix>.tsv analyse.R's columns plus the per-gene predictions
  causal<suffix>.tsv      ID gene beta MAC score_informative per causal variant
  sim_meta<suffix>.json   the knobs, the realised h2 and the prediction table
usage: make_power_pheno.py --dir <fixture> [--scores <dir>] [--suffix _p2] ...
"""
import argparse, collections, gzip, json, os, sys, time
import numpy as np
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True)
ap.add_argument("--scores", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "..", "..", "work", "flexrv_study", "scores"))
ap.add_argument("--suffix", default="_p2")
ap.add_argument("--seed", type=int, default=20260914)
ap.add_argument("--n-causal", type=int, default=60)
ap.add_argument("--ladder", default="10,80", help="log-uniform NCP of the causal-only burden")
ap.add_argument("--maxmaf", type=float, default=0.001, help="the study cell's max MAF")
ap.add_argument("--min-totmac", type=int, default=30, help="eligibility: visible-cell minor alleles")
ap.add_argument("--k", type=float, default=10.0, help="logistic slope of P(causal) in the latent pct")
ap.add_argument("--pmax", type=float, default=0.9, help="P(causal) at the top of the latent")
ap.add_argument("--prevalence", type=float, default=0.10)
ap.add_argument("--dry-run", action="store_true", help="print the prediction table only")
a = ap.parse_args()
t0 = time.time()
lo, hi = (float(x) for x in a.ladder.split(","))
assert 0 < lo < hi, a.ladder


def say(*m):
    print(f"[{time.time() - t0:6.1f}s]", *m, flush=True)


# ------------------------------------------------------------ inputs
V = []   # (row, ID, gene, anno, AF, MAF, AC, a_lat) over the coding set, bim order
with gzip.open(os.path.join(a.dir, "variants.tsv.gz"), "rt") as fh:
    h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
    for row, line in enumerate(fh):
        f = line.rstrip("\n").split("\t")
        if f[ix["gene"]] == ".":
            continue
        V.append((row, f[ix["ID"]], f[ix["gene"]], f[ix["anno"]], float(f[ix["AF"]]),
                  float(f[ix["MAF"]]), int(float(f[ix["AC"]])), float(f[ix["a_lat"]])))
bim_ids = [l.split()[1] for l in open(os.path.join(a.dir, "chr20_sim.bim"))]
for row, vid, *_ in V[::997]:
    assert bim_ids[row] == vid, f"bim/tsv order differs at row {row}"
meta = json.load(open(os.path.join(a.dir, "sim_meta.json")))
N = meta["N"]; NG = meta["n_genes"]
say(f"{len(V):,} coding variants, N={N}, {NG} genes, seed {meta['seed']} (fixture) / {a.seed} (this)")


def read_scores(arm):
    """ID -> score from the arm's group file (var line then anno then score:SIM)."""
    out = {}
    lines = open(os.path.join(a.scores, f"groupFile_{arm}.txt")).read().splitlines()
    i = 0
    while i < len(lines):
        v = lines[i].split(); assert v[1] == "var", lines[i][:40]
        j = i + 1
        while lines[j].split()[1] != "score:SIM":
            j += 1
        s = lines[j].split()[2:]; assert len(s) == len(v) - 2, v[0]
        out.update(zip(v[2:], (float(x) for x in s)))
        i = j + 1
    return out


scores = {arm: read_scores(arm) for arm in ("informative", "realistic", "uninformative")}
say("scores read: " + ", ".join(f"{k} {len(v):,}" for k, v in scores.items()))


def pct(x):
    r = np.argsort(np.argsort(x, kind="stable"), kind="stable") + 1.0
    return r / len(x)


# ------------------------------------------------------------ the visible cell
genes = collections.defaultdict(list)
for rec in V:
    if rec[3] in ("high", "medium") and rec[5] <= a.maxmaf:
        genes[rec[2]].append(rec)
totmac = {g: sum(min(r[6], 2 * N - r[6]) for r in v) for g, v in genes.items()}
elig = sorted(g for g in genes if totmac[g] >= a.min_totmac)
rng = np.random.default_rng(a.seed)
cg = sorted(rng.choice(elig, min(a.n_causal, len(elig)), replace=False))
target = dict(zip(cg, np.exp(rng.uniform(np.log(lo), np.log(hi), len(cg)))))
say(f"visible cell: {sum(len(v) for v in genes.values()):,} variants in {len(genes)} genes; "
    f"{len(elig)} eligible (>= {a.min_totmac} minor alleles); {len(cg)} causal genes drawn")

# ------------------------------------------------------------ causal variants + prediction
thr = stats.chi2.isf(0.05 / NG, 1)
EXP = (0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20)
causal_rows = []      # (row, ID, gene, beta, MAC, s_inf)
pred = collections.defaultdict(list)
gene_pred = {}
sign = {}
for g in cg:
    v = genes[g]
    maf = np.array([r[5] for r in v]); al = np.array([r[7] for r in v])
    hi_cls = np.array([r[3] == "high" for r in v]); nh = ~hi_cls
    z = np.full(len(v), 0.6)
    if nh.any():
        z[nh] = pct(al[nh]) - 0.5
    p = a.pmax / (1.0 + np.exp(-a.k * z))
    causal = rng.random(len(v)) < p
    if not causal.any():
        causal[np.argmax(p)] = True
    w2 = 2 * maf * (1 - maf)
    Sc = w2[causal].sum()
    b = np.sqrt(target[g] / (N * Sc))            # HWE prediction; delivered exactly below
    sign[g] = 1.0 if rng.random() < 0.5 else -1.0
    ids = [r[1] for r in v]
    sc = {arm: np.array([scores[arm][i] for i in ids]) for arm in scores}

    def ncp(w):
        return N * ((w * w2 * b * causal).sum()) ** 2 / ((w ** 2 * w2).sum())

    gp = {"unw": ncp(np.ones(len(v))), "orc": ncp(causal.astype(float))}
    for arm, s in sc.items():
        gp[f"{arm}_s1"] = ncp(s)
        gp[f"{arm}_s3"] = ncp(s ** 3)
        gp[f"{arm}_best"] = max(ncp(s ** e) for e in EXP)
    for k_, x in gp.items():
        pred[k_].append(x)
    gene_pred[g] = dict(gp, n_vis=len(v), n_causal=int(causal.sum()),
                        cmac=int(sum(min(r[6], 2 * N - r[6]) for r, c in zip(v, causal) if c)),
                        cor_inf=float(np.corrcoef(sc["informative"], causal)[0, 1]) if causal.sum() < len(v) else float("nan"))
    for r, c, s_inf in zip(v, causal, sc["informative"]):
        if c:
            causal_rows.append((r[0], r[1], g, sign[g] * b, min(r[6], 2 * N - r[6]), s_inf))


def edet(x, f=1.0):
    return sum(stats.ncx2.sf(thr, 1, v * f) for v in x)


print(f"\nPREDICTION (before generation): {len(cg)} causal genes, ladder [{lo}, {hi}] on the causal-only "
      f"burden, Bonferroni threshold chi2 = {thr:.2f} (alpha 0.05/{NG})")
print(f"  causal variants/gene median {np.median([gene_pred[g]['n_causal'] for g in cg]):.0f}, "
      f"visible/gene median {np.median([gene_pred[g]['n_vis'] for g in cg]):.0f}, "
      f"causal minor alleles/gene median {np.median([gene_pred[g]['cmac'] for g in cg]):.0f} "
      f"(q10 {np.percentile([gene_pred[g]['cmac'] for g in cg], 10):.0f}), "
      f"mean cor(informative score, causal) {np.nanmean([gene_pred[g]['cor_inf'] for g in cg]):.2f}")
print(f"  {'burden':20s} NCP median   E[detections] (quant = binary by construction)")
for k_ in ("unw", "orc", "informative_s1", "informative_s3", "informative_best", "realistic_s1",
           "realistic_s3", "realistic_best", "uninformative_s1", "uninformative_s3", "uninformative_best"):
    print(f"  {k_:20s} {np.median(pred[k_]):8.1f}   {edet(pred[k_]):6.1f}")
print("  '_best' is the max over exponents 0.2..20 of the score, an UPPER bound on what a Cauchy\n"
      "  combination over the transform sets can return; the study's uninformative arm is expected\n"
      "  between unw and uninformative_best, and the informative arm between informative_s1 and\n"
      "  informative_best.\n")
if a.dry_run:
    sys.exit(0)

# ------------------------------------------------------------ genotypes of the causal set
bed = os.path.join(a.dir, "chr20_sim.bed")
nb = (N + 3) // 4
unpack = np.array([2, np.nan, 1, 0])           # plink1 code -> ALT dosage (02_build_dataset.py)
g_cod = np.zeros(N)
burden = collections.defaultdict(lambda: np.zeros(N))
n_missing = 0
af_by_row = {r[0]: r[4] for r in V}
with open(bed, "rb") as fh:
    for row, vid, g, beta, mac, _ in sorted(causal_rows):
        fh.seek(3 + row * nb)
        raw = np.frombuffer(fh.read(nb), np.uint8)
        codes = np.empty((nb, 4), np.uint8)
        for s in range(4):
            codes[:, s] = (raw >> (2 * s)) & 3
        d = unpack[codes.ravel()[:N]]
        m = np.isnan(d); n_missing += int(m.sum()); d[m] = 0.0
        if af_by_row[row] > 0.5:
            d = 2.0 - d                          # the test folds to the minor allele
        burden[g] += beta * d
say(f"{len(causal_rows)} causal variants read from the .bed ({n_missing} missing calls set to 0)")
g_cod = sum(burden.values())

# ------------------------------------------------------------ scale + traits
# The ladder is delivered on the REALISED per-gene burden variance (relatedness
# and LD move it off the HWE value used above): per gene, scale beta so that
# N * var(burden_g) equals the drawn target exactly, then the total h2 is the
# sum, an outcome.
for g in cg:
    v = burden[g].var()
    f = np.sqrt(target[g] / (N * v)) if v > 0 else 0.0
    burden[g] *= f
    sign[g] = (sign[g], f)
g_cod = sum(burden.values())
H2 = float(g_cod.var())
tq = stats.norm.isf(a.prevalence)
BIN_NCP_FRAC = float(stats.norm.pdf(tq) ** 2 / (a.prevalence * (1 - a.prevalence)))
say(f"realised chr20 h2 (quant) = {H2:.4f}; K = {a.prevalence}, linearised NCP fraction {BIN_NCP_FRAC:.4f} "
    f"(reference only: the binary liability is solved per gene on the exact threshold model below)")

ph = open(os.path.join(a.dir, "pheno_pcs.txt")).read().splitlines()
hdr = ph[0].split(); iy = hdr.index("y_null")
y_null = np.array([float(l.split()[iy]) for l in ph[1:]])
assert len(y_null) == N
fam_iids = [l.split()[1] for l in open(os.path.join(a.dir, "chr20_sim.fam"))]
assert fam_iids == [l.split()[hdr.index("IID")] for l in ph[1:]], "pheno rows are not in fam order"
gs = g_cod / np.sqrt(H2)                               # unit variance
yq = np.sqrt(H2) * gs + np.sqrt(1 - H2) * y_null

# ---- the binary arm: one scale c_g PER GENE on the exact liability-threshold
# model (2026-09-14 evening). The first version scaled every gene by the same
# sqrt(h2 / BIN_NCP_FRAC), the linearised constant, and the sweep found the
# binary arm at 0.22 of its designed non-centrality above ~0.6 sd per allele
# (handoff-2026-09-14-flexrv-p2.md, 4b): a carrier shifted by 1.7 sd is a case
# with probability ~0.5, not K plus a linear increment. The score test's
# non-centrality for burden B on a 0/1 trait is N Cov(y, B)^2 / (K(1-K) Var(B))
# with P(case | B = b) = Phi((c b - t) / sd_rest) under the threshold model,
# which saturates in c; c_g is solved so that it equals the gene's quantitative
# NCP, and the fixed point in the residual sd is iterated. A gene whose target
# is unreachable (every carrier a case is the ceiling, ~ (1-K)/K * MAC) is set
# to 0.9 of its ceiling and flagged in truth_genes<suffix>.tsv. Every binary
# effect is RISK-INCREASING: the sign is immaterial to the quantitative test,
# and on a 0/1 trait a protective rare burden has a ceiling K/(1-K) times
# lower (carriers can at best never be cases), which no ladder could reach.
from scipy.optimize import brentq
Bpos = {g: burden[g] * sign[g][0] for g in cg}     # non-negative burdens
K = a.prevalence
def ncp_exact(c, B, t, sd_rest):
    pc = stats.norm.cdf((c * B - t) / sd_rest)
    cov = np.mean(B * pc) - np.mean(B) * np.mean(pc)
    return N * cov ** 2 / (K * (1 - K) * B.var())
c_bin = {g: 1.0 for g in cg}; unreachable = set(); ncp_bin = {}
for it in range(3):
    gb = sum(c_bin[g] * Bpos[g] for g in cg)
    Vb = float(gb.var())
    liab_b = gb + np.sqrt(max(1e-9, 1 - Vb)) * y_null
    t = float(np.quantile(liab_b, 1 - K))
    for g in cg:
        B = Bpos[g]
        sd_rest = np.sqrt(max(1e-9, liab_b.var() - c_bin[g] ** 2 * B.var()))
        cmax = 12.0 / B.max()
        ceiling = ncp_exact(cmax, B, t, sd_rest)
        tgt = target[g]
        if ceiling <= tgt:
            unreachable.add(g); tgt = 0.9 * ceiling
        c_bin[g] = brentq(lambda c: ncp_exact(c, B, t, sd_rest) - tgt, 0.0, cmax)
        ncp_bin[g] = ncp_exact(c_bin[g], B, t, sd_rest)
gb = sum(c_bin[g] * Bpos[g] for g in cg)
h2b = float(gb.var())
liab_b = gb + np.sqrt(1 - h2b) * y_null
yb = (liab_b > np.quantile(liab_b, 1 - K)).astype(int)
# two ways: the realised marginal score statistic per gene, Z = U / sqrt(Var0 U)
# with U = sum (y - ybar)(B - Bbar) and Var0 the NULL variance, against
# sqrt(NCP). E[Z] = sqrt(NCP) on both arms; E[Z^2] - 1 is NOT the NCP on the
# binary arm because the alternative variance of U exceeds the null one (a
# carrier with P(case) ~ 0.5 contributes 0.25, not K(1-K)), so the check is
# on the mean of Z, whose noise over 60 genes is ~0.13.
def zreal(y, Bd):
    return np.array([np.sum((y - y.mean()) * (Bd[g] - Bd[g].mean()))
                     / np.sqrt(y.var() * np.sum((Bd[g] - Bd[g].mean()) ** 2)) for g in cg])
zb, zq = zreal(yb, Bpos), zreal(yq, burden)      # the quant arm's burdens are signed
sq_b = np.sqrt([ncp_bin[g] for g in cg]); sq_q = np.sqrt([target[g] for g in cg])
say(f"binary arm solved per gene on the exact threshold model: liability h2 {h2b:.4f}, "
    f"{len(unreachable)} of {len(cg)} genes unreachable (set to 0.9 x ceiling); "
    f"delivered NCP_binary median {np.median([ncp_bin[g] for g in cg]):.1f} vs NCP_quant median "
    f"{np.median([target[g] for g in cg]):.1f}; realised mean Z / mean sqrt(NCP): binary "
    f"{zb.mean() / sq_b.mean():.3f}, quant {zq.mean() / sq_q.mean():.3f} (expect 1 +- 0.03)")
ynb = (y_null > np.quantile(y_null, 1 - K)).astype(int)
say(f"y_quant{a.suffix}: var {yq.var():.3f}; y_binary{a.suffix}: {yb.sum()} cases "
    f"({yb.mean() * 100:.2f}%); y_null_bin{a.suffix}: {ynb.sum()} cases")

# ------------------------------------------------------------ outputs
sfx = a.suffix
with open(os.path.join(a.dir, f"pheno_pcs{sfx}.txt"), "w") as fh:
    fh.write(" ".join(hdr + [f"y_quant{sfx}", f"y_binary{sfx}", f"y_null_bin{sfx}"]) + "\n")
    for i, l in enumerate(ph[1:]):
        fh.write(f"{l} {yq[i]:.6f} {yb[i]} {ynb[i]}\n")

with open(os.path.join(a.dir, f"causal{sfx}.tsv"), "w") as fh:
    fh.write("ID\tgene\tbeta\tMAC\tscore_informative\n")
    for row, vid, g, beta, mac, s_inf in sorted(causal_rows):
        fh.write(f"{vid}\t{g}\t{beta * sign[g][1]:.6g}\t{mac}\t{s_inf:.4f}\n")

old = {l.split("\t")[0]: l.split("\t") for l in open(os.path.join(a.dir, "truth_genes.tsv")).read().splitlines()[1:]}
with open(os.path.join(a.dir, f"truth_genes{sfx}.tsv"), "w") as fh:
    fh.write("gene\tcausal\tn_var\tn_high\tn_medium\tn_low\tn_causal_var\th2\tNCP_target\tNCP_quant\t"
             "NCP_binary\tdirectional\tdir_index\tlink\tn_visible\tcausal_MAC\tpred_unw\tpred_orc\t"
             "pred_informative_s1\tpred_informative_best\tpred_realistic_s1\tpred_realistic_best\t"
             "pred_uninformative_s1\tpred_uninformative_best\tc_bin\tbinary_unreachable\n")
    for g in sorted(old):
        o = old[g]
        if g in gene_pred:
            gp = gene_pred[g]; nq = N * burden[g].var()
            fh.write(f"{g}\t1\t{o[2]}\t{o[3]}\t{o[4]}\t{o[5]}\t{gp['n_causal']}\t{burden[g].var():.4e}\t"
                     f"{target[g]:.2f}\t{nq:.2f}\t{ncp_bin[g]:.2f}\t1\t1.000\tgraded\t{gp['n_vis']}\t{gp['cmac']}\t"
                     f"{gp['unw']:.2f}\t{gp['orc']:.2f}\t{gp['informative_s1']:.2f}\t{gp['informative_best']:.2f}\t"
                     f"{gp['realistic_s1']:.2f}\t{gp['realistic_best']:.2f}\t{gp['uninformative_s1']:.2f}\t"
                     f"{gp['uninformative_best']:.2f}\t{c_bin[g]:.4f}\t{int(g in unreachable)}\n")
        else:
            fh.write(f"{g}\t0\t{o[2]}\t{o[3]}\t{o[4]}\t{o[5]}\t0\t0.0000e+00\t0.00\t0.00\t0.00\t0\t0.000\t.\t"
                     f"{len(genes.get(g, []))}\t0" + "\t0.00" * 8 + "\t0.0000\t0\n")

json.dump({"suffix": sfx, "seed": a.seed, "n_causal_genes": len(cg), "ladder": a.ladder,
           "maxmaf": a.maxmaf, "min_totmac": a.min_totmac, "k": a.k, "pmax": a.pmax,
           "prevalence": a.prevalence, "bin_ncp_frac": BIN_NCP_FRAC, "h2_realised": H2,
           "h2_binary_liability": h2b, "binary_model": "exact threshold, c_g per gene",
           "binary_unreachable_genes": sorted(unreachable),
           "n_causal_variants": len(causal_rows), "n_missing_calls": n_missing,
           "background": "y_null of pheno_pcs.txt, scaled by sqrt(1 - h2)",
           "expected_detections": {k_: edet(v) for k_, v in pred.items()},
           "bonferroni_chi2": thr}, open(os.path.join(a.dir, f"sim_meta{sfx}.json"), "w"), indent=2)
say("done")
