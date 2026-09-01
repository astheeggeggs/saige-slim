#!/usr/bin/env python
"""Turn the chr20 tree sequence into a PLINK dataset + traits + annotations.

Pipeline
  1. Take the EUR individuals only -- the AFR samples exist so the published
     two-population model is simulated as published; they are dropped here.
  2. Build 500 nuclear families by gene-drop from 1000 of the CEU founders, so
     the cohort has real relatedness for the sparse GRM to find (a pure
     coalescent sample has none: every off-diagonal would be ~0 and step 0
     would return a diagonal matrix).
  3. Stream every biallelic site, write .bed/.bim/.fam directly.
  4. Lay ~550 genes with exon structure over chr20; annotate coding variants
     into low/medium/high deleteriousness classes that are FUZZY readouts of a
     latent score which itself depends on rarity -- so high/medium are rarer on
     average, and annotation predicts effect only imperfectly.
  5. Spike-and-slab effects with MAF-dependent slab variance (alpha).
  6. Continuous + binary (liability threshold) traits.

Effect-size model
  Latent deleteriousness for coding variant j:
      d_j = gamma * z(-log10 MAF_j) + N(0,1)
  Two noisy readouts of the same d_j:
      annotation  a_j = d_j + N(0, sigma_a^2)   -> quantile-cut to low/med/high
      true class  e_j = d_j + N(0, sigma_e^2)   -> quantile-cut to null/bi/dir
  so cor(a,e) < 1 ("fuzzy but correlated") and both track rarity.
  Slab sd scales as [2p(1-p)]^(alpha/2), i.e. beta_j ~ N(0, s^2 [2p(1-p)]^alpha),
  the standard alpha-model. With alpha=-0.3 rarer variants get larger effects
  but per-variant h2 still scales as [2p(1-p)]^(1+alpha).
  Only a random 10% of genes are causal; the other 90% are exactly null and are
  the built-in negative control for the gene-based tests.
"""
import argparse
import gzip
import json
import math
import os
import statistics
import time

import numpy as np
import stdpopsim
import tskit

ap = argparse.ArgumentParser()
ap.add_argument("--trees", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--prefix", default="chr20_sim")
ap.add_argument("--eur-pop", default="EUR")
ap.add_argument("--seed", type=int, default=7)
ap.add_argument("--n-singletons", type=int, default=8000)
ap.add_argument("--n-families", type=int, default=500)
ap.add_argument("--n-genes", type=int, default=550)
ap.add_argument("--frac-causal-genes", type=float, default=0.10)
ap.add_argument("--alpha", type=float, default=-0.3, help="MAF-dependence of the slab")
ap.add_argument("--h2", type=float, default=0.01, help="chr20 h2 of the primary traits")
ap.add_argument("--h2-pc", type=float, default=0.20, help="chr20 h2 of the positive controls")
ap.add_argument("--h2-poly", type=float, default=0.35,
                help="background polygenic h2 from the UNSIMULATED rest of the genome, "
                     "transmitted through the pedigree")
ap.add_argument("--frac-h2-coding", type=float, default=0.6)
ap.add_argument("--prevalence", type=float, default=0.10)
ap.add_argument("--var-covar", type=float, default=0.10, help="variance explained by age+sex")
ap.add_argument("--ncp-ladder", default="",
                help='"lo,hi": give each causal gene a target non-centrality drawn '
                     "LOG-UNIFORMLY on [lo,hi] and scale that gene's effects to hit "
                     "it, instead of scaling all coding effects to a single --h2. "
                     "Power then comes out as a graded curve and the assertion "
                     "becomes Spearman(NCP, -log10 p) over genes rather than a hit "
                     "count, which is the only thing that still discriminates once "
                     "N is large enough to saturate a fixed h2. --h2 then sets the "
                     "NON-coding component only; total h2 is an outcome, reported.")
ap.add_argument("--h2-binary-mult", default="1.0",
                help="multiply the genetic component by this for the BINARY arm's "
                     "liability. A binary trait at prevalence K buys far less "
                     "non-centrality than a quantitative one at the same h2: "
                     "NCP_bin/NCP_quant = phi(t)^2/(K(1-K)), which is 0.0717 at "
                     "K=0.01, so the two arms need separately scaled heritability "
                     "or the binary ladder sits ~14x below the quantitative one. "
                     "\"auto\" solves for the multiplier that equalises the two "
                     "arms at the current --prevalence, so it cannot go stale when "
                     "the prevalence changes; a number sets it directly.")
ap.add_argument("--thin-intergenic", type=float, default=1.0,
                help="retention probability for RARE INTERGENIC variants (1.0 = keep "
                     "everything, the default). Thinning must never be uniform: the "
                     "GRM (--minMAFforGRM) and the LD-pruned variance-ratio set are "
                     "drawn from the common-variant pool, and the group tests from the "
                     "genic one, so a uniform thin destroys exactly the pools the "
                     "fixture exists to exercise. Genic variants and variants at or "
                     "above --thin-maf are always kept.")
ap.add_argument("--thin-maf", type=float, default=0.01,
                help="MAF at or above which a variant is never thinned")
a = ap.parse_args()

rng = np.random.default_rng(a.seed)
# Separate stream: with --thin-intergenic at its 1.0 default no draw is taken
# and the main stream is bit-identical to an unthinned run, so turning thinning
# on changes which variants exist WITHOUT also reshuffling every effect size.
rng_thin = np.random.default_rng(a.seed + 90001)
t0 = time.time()


def say(*m):
    print(f"[{time.time() - t0:7.1f}s]", *m, flush=True)


# A binary trait at prevalence K buys phi(t)^2 / (K(1-K)) of the non-centrality
# a quantitative trait gets from the same liability h2 (t = qnorm(1-K)). At
# K=0.01 that is 0.0718, i.e. the binary arm is 13.9x weaker unless its
# heritability is scaled up to compensate -- which is what --h2-binary-mult does.
_t = statistics.NormalDist().inv_cdf(1 - a.prevalence)
BIN_NCP_FRAC = (math.exp(-_t * _t / 2) / math.sqrt(2 * math.pi)) ** 2 / (
    a.prevalence * (1 - a.prevalence))
if str(a.h2_binary_mult).lower() == "auto":
    a.h2_binary_mult = 1.0 / BIN_NCP_FRAC
else:
    a.h2_binary_mult = float(a.h2_binary_mult)
say(f"prevalence {a.prevalence}: binary/quantitative NCP ratio {BIN_NCP_FRAC:.4f} "
    f"at equal liability h2; --h2-binary-mult={a.h2_binary_mult:.4g} "
    f"-> net {a.h2_binary_mult * BIN_NCP_FRAC:.3f}x")

# ---------------------------------------------------------------- 1. samples
ts = tskit.load(a.trees)
pop = {p.metadata["name"]: p.id for p in ts.populations()}
node_pop = ts.tables.nodes.population
ceu = [ind.id for ind in ts.individuals() if node_pop[ind.nodes[0]] == pop[a.eur_pop]]
n_founder = a.n_singletons + 2 * a.n_families
assert len(ceu) >= n_founder, (len(ceu), n_founder)
founders = np.array(ceu[:n_founder])
founder_nodes = np.concatenate([ts.individual(i).nodes for i in founders]).astype(np.int32)
N = a.n_singletons + 4 * a.n_families
n_child = 2 * a.n_families
say(f"{a.eur_pop} individuals {len(ceu)}; founders {n_founder}; output N {N} "
    f"({a.n_singletons} unrelated + {a.n_families} families of 4)")

# --------------------------------------------- 2. meiosis crossover positions
gmap = stdpopsim.get_species("HomSap").get_contig(
    "chr20", genetic_map="HapMapII_GRCh38").recombination_map
gpos = np.asarray(gmap.position)
grate = np.nan_to_num(np.asarray(gmap.rate), nan=0.0)
gcum = np.concatenate([[0.0], np.cumsum(grate * np.diff(gpos))])
total_morgans = gcum[-1]
say(f"chr20 genetic length {total_morgans * 100:.1f} cM")

# transmissions 0..n_child-1 = paternal, n_child..2*n_child-1 = maternal
n_trans = 2 * n_child
n_xo = rng.poisson(total_morgans, n_trans)
ev_pos, ev_trans = [], []
for t in range(n_trans):
    if n_xo[t]:
        u = np.sort(rng.uniform(0, total_morgans, n_xo[t]))
        ev_pos.append(np.interp(u, gcum, gpos))
        ev_trans.append(np.full(n_xo[t], t))
ev_pos = np.concatenate(ev_pos) if ev_pos else np.zeros(0)
ev_trans = np.concatenate(ev_trans).astype(np.int32) if len(ev_trans) else np.zeros(0, np.int32)
o = np.argsort(ev_pos, kind="stable")
ev_pos, ev_trans = ev_pos[o], ev_trans[o]
cur_hap = rng.integers(0, 2, n_trans).astype(np.int8)
say(f"{n_trans} transmissions, {len(ev_pos)} crossovers "
    f"({len(ev_pos) / n_trans:.2f} per meiosis)")

# founder-node index of each child's father / mother (base of the 2-node pair)
fam_of_child = np.repeat(np.arange(a.n_families), 2)
fa_base = 2 * (a.n_singletons + 2 * fam_of_child)
mo_base = 2 * (a.n_singletons + 2 * fam_of_child + 1)
# output-row indices
row_par = (a.n_singletons + 4 * np.arange(a.n_families)[:, None] + np.array([0, 1])).ravel()
row_chi = (a.n_singletons + 4 * np.arange(a.n_families)[:, None] + np.array([2, 3])).ravel()

# ------------------------------------------------------------------ 3. genes
# Place genes over the region that actually carries variation, not over the
# nominal sequence length: with --right the ts keeps chr20's full length but
# only the truncated prefix is polymorphic, so a uniform draw over L would put
# almost every gene in an empty region.
sp = ts.tables.sites.position
lo, hi = (float(sp.min()), float(sp.max())) if len(sp) else (5e5, ts.sequence_length)
L = ts.sequence_length
starts = np.sort(rng.uniform(lo, hi, a.n_genes))
spans = np.clip(rng.lognormal(np.log(30000), 0.9, a.n_genes), 8000, 3e5)
spans = np.minimum(spans, np.diff(np.append(starts, hi)) * 0.9)
cds = np.clip(rng.lognormal(np.log(2000), 0.8, a.n_genes), 300, 15000)
ex_s, ex_e, ex_g = [], [], []
for g in range(a.n_genes):
    ne = max(1, int(round(cds[g] / 160)))
    elen = cds[g] / ne
    slots = np.linspace(starts[g], starts[g] + spans[g] - elen, ne)
    jit = rng.uniform(0, max(1.0, (spans[g] - cds[g]) / max(ne, 1) * 0.5), ne)
    s = np.sort(slots + jit)
    ex_s.append(s)
    ex_e.append(s + elen)
    ex_g.append(np.full(ne, g))
ex_s = np.concatenate(ex_s)
ex_e = np.concatenate(ex_e)
ex_g = np.concatenate(ex_g).astype(np.int32)
oo = np.argsort(ex_s)
ex_s, ex_e, ex_g = ex_s[oo], ex_e[oo], ex_g[oo]
say(f"{a.n_genes} genes, {len(ex_s)} exons, {ex_e.sum() - ex_s.sum():,.0f} coding bp "
    f"({(ex_e - ex_s).sum() / L * 100:.2f}% of chr20)")

# -------------------------------------------------- 4. stream sites -> .bed
CODE = np.array([3, 2, 0], np.uint8)  # dosage 0,1,2 -> plink1 0b11,0b10,0b00
bed = open(f"{a.outdir}/{a.prefix}.bed", "wb")
bed.write(bytes([0x6C, 0x1B, 0x01]))
bim = open(f"{a.outdir}/{a.prefix}.bim", "w")

buf, BUFMAX = [], 4000
pos_l, ac_l, gene_l = [], [], []
cod_idx = []
# Coding dosages go straight to disk as raw int8 and are read back as a memmap.
# Held as a list of N-vectors they cost n_coding * N bytes resident, which
# np.stack then briefly DOUBLES -- 2.5 GB each at N=50,000.
cod_path = f"{a.outdir}/_coding_dosage.i8"
cod_fh = open(cod_path, "wb")
ev_ptr = 0
n_seen = n_kept = n_multi = n_thinned = 0
thin = a.thin_intergenic

for var in ts.variants(samples=founder_nodes, copy=False):
    n_seen += 1
    if len(var.alleles) != 2 or var.alleles[1] is None:
        n_multi += 1
        continue
    p = var.site.position
    while ev_ptr < len(ev_pos) and ev_pos[ev_ptr] <= p:
        cur_hap[ev_trans[ev_ptr]] ^= 1
        ev_ptr += 1
    g = var.genotypes
    dose = np.empty(N, np.int8)
    dose[:a.n_singletons] = (g[:2 * a.n_singletons].reshape(-1, 2).sum(1))
    dose[row_par] = g[2 * a.n_singletons:].reshape(-1, 2).sum(1)
    dose[row_chi] = (g[fa_base + cur_hap[:n_child]].astype(np.int8)
                     + g[mo_base + cur_hap[n_child:]].astype(np.int8))
    ac = int(dose.sum())
    if ac == 0 or ac == 2 * N:
        continue
    k = np.searchsorted(ex_s, p, "right") - 1
    gi = int(ex_g[k]) if (k >= 0 and p < ex_e[k]) else -1
    if thin < 1.0 and gi < 0:
        af_i = ac / (2 * N)
        if min(af_i, 1 - af_i) < a.thin_maf and rng_thin.random() >= thin:
            n_thinned += 1
            continue
    codes = CODE[dose].reshape(-1, 4)
    buf.append(codes[:, 0] | (codes[:, 1] << 2) | (codes[:, 2] << 4) | (codes[:, 3] << 6))
    if len(buf) >= BUFMAX:
        bed.write(np.concatenate(buf).tobytes())
        buf = []
    ip = int(p)
    ref, alt = var.alleles[0], var.alleles[1]
    bim.write(f"20\t20:{ip}:{ref}:{alt}\t0\t{ip}\t{alt}\t{ref}\n")
    pos_l.append(ip)
    ac_l.append(ac)
    gene_l.append(gi)
    if gi >= 0:
        cod_fh.write(dose.tobytes())
        cod_idx.append(n_kept)
    n_kept += 1

if buf:
    bed.write(np.concatenate(buf).tobytes())
bed.close()
bim.close()
cod_fh.close()
say(f"sites seen {n_seen:,} (dropped {n_multi:,} non-biallelic); "
    f"kept {n_kept:,} polymorphic in EUR; coding {len(cod_idx):,}")
if thin < 1.0:
    say(f"thinned {n_thinned:,} intergenic variants with MAF < {a.thin_maf} "
        f"(retention {thin}); genic and MAF >= {a.thin_maf} untouched")

pos = np.array(pos_l)
ac = np.array(ac_l, np.int64)
gene = np.array(gene_l, np.int32)
af = ac / (2 * N)
maf = np.minimum(af, 1 - af)
del pos_l, ac_l, gene_l

# ---------------------------------------------------------------- .fam file
iid = np.array([f"S{i + 1:06d}" for i in range(N)])
fid = np.array(["UNREL"] * N, dtype=object)
pat = np.array(["0"] * N, dtype=object)
mat = np.array(["0"] * N, dtype=object)
for f in range(a.n_families):
    b = a.n_singletons + 4 * f
    for k in range(4):
        fid[b + k] = f"FAM{f + 1:04d}"
    pat[b + 2] = pat[b + 3] = iid[b]
    mat[b + 2] = mat[b + 3] = iid[b + 1]
sex = rng.integers(1, 3, N)          # 1=male 2=female (plink coding)
sex[row_par[0::2]] = 1               # fathers male
sex[row_par[1::2]] = 2               # mothers female
with open(f"{a.outdir}/{a.prefix}.fam", "w") as fh:
    for i in range(N):
        fh.write(f"{fid[i]} {iid[i]} {pat[i]} {mat[i]} {sex[i]} -9\n")

# Pedigree truth: every nonzero off-diagonal of the relatedness matrix that
# step 0 is trying to estimate. Founders are drawn from the coalescent and the
# two parents of a family are unrelated, so first-degree pairs within a family
# are the complete list. With this on disk, step 0 becomes a CONTROL -- compare
# its .mtx against this -- instead of an unchecked dependency of everything
# downstream. "grm" is the SAIGE scale (2 * kinship: 1.0 self, 0.5 first-degree).
with open(f"{a.outdir}/{a.prefix}.pedigree_truth.txt", "w") as fh:
    fh.write("IID1\tIID2\tkinship\tgrm\trelation\n")
    for i in range(N):
        fh.write(f"{iid[i]}\t{iid[i]}\t0.5\t1.0\tself\n")
    n_pair = 0
    for f in range(a.n_families):
        b = a.n_singletons + 4 * f
        fa, mo, c1, c2 = iid[b], iid[b + 1], iid[b + 2], iid[b + 3]
        for x, y, rel in ((fa, c1, "parent-offspring"), (fa, c2, "parent-offspring"),
                          (mo, c1, "parent-offspring"), (mo, c2, "parent-offspring"),
                          (c1, c2, "full-sib")):
            fh.write(f"{x}\t{y}\t0.25\t0.5\t{rel}\n")
            n_pair += 1
say(f"pedigree truth: {N} self + {n_pair:,} first-degree pairs")

np.save(f"{a.outdir}/_pos.npy", pos)
np.save(f"{a.outdir}/_ac.npy", ac)
np.save(f"{a.outdir}/_gene.npy", gene)
bim_ids = np.array([l.split("\t")[1] for l in
                    open(f"{a.outdir}/{a.prefix}.bim")], dtype=object)
assert len(bim_ids) == n_kept
say("wrote bed/bim/fam")

# ------------------------------------------- 5. annotations + effect classes
cod_idx = np.array(cod_idx)
D = (np.memmap(cod_path, np.int8, "r", shape=(len(cod_idx), N))
     if len(cod_idx) else np.zeros((0, N), np.int8))
cmaf = maf[cod_idx]
u = -np.log10(cmaf)
z = (u - u.mean()) / u.std()
GAMMA, SIG_A, SIG_E = 0.9, 0.8, 0.8
d = GAMMA * z + rng.normal(size=len(z))
a_lat = d + rng.normal(0, SIG_A, len(z))
e_lat = d + rng.normal(0, SIG_E, len(z))

P_HIGH, P_MED = 0.05, 0.30
def cut(x, p_hi, p_mid):
    q = np.quantile(x, [1 - p_hi - p_mid, 1 - p_hi])
    return np.digitize(x, q)          # 0 = low, 1 = medium, 2 = high
anno = cut(a_lat, P_HIGH, P_MED)
eclass = cut(e_lat, P_HIGH, P_MED)
ANN = np.array(["low", "medium", "high"])
say("annotation x MAF: " + ", ".join(
    f"{ANN[c]} n={int((anno == c).sum())} medianMAF={np.median(cmaf[anno == c]):.2e}"
    for c in (2, 1, 0)))
say(f"cor(anno, true effect class) = {np.corrcoef(a_lat, e_lat)[0, 1]:.3f} (latent), "
    f"{np.corrcoef(anno, eclass)[0, 1]:.3f} (after cutting)")

# ------------------------------------------------------ 6. spike-and-slab
gcode = gene[cod_idx]
n_var_per_gene = np.bincount(gcode[gcode >= 0], minlength=a.n_genes)
# The group file skips genes with < 2 coding variants, so such a gene cannot be
# tested at all; making it "causal" only pollutes the positive-control set.
eligible = np.where(n_var_per_gene >= 2)[0]
n_causal_genes = int(round(a.frac_causal_genes * a.n_genes))
if n_causal_genes > len(eligible):
    print(f"WARNING: only {len(eligible)} genes have >= 2 coding variants; "
          f"lowering causal genes {n_causal_genes} -> {len(eligible)}", flush=True)
    n_causal_genes = len(eligible)
causal_genes = np.sort(rng.choice(eligible, n_causal_genes, replace=False))
is_cg = np.isin(gcode, causal_genes)
scale = (2 * cmaf * (1 - cmaf)) ** (a.alpha / 2)      # slab sd multiplier

PI = np.array([0.02, 0.35, 0.55])                     # P(causal | true class)
SD = np.array([0.25, 1.00, 1.20])                     # slab sd  | true class
beta_c = np.zeros(len(cod_idx))
causal = (rng.random(len(cod_idx)) < PI[eclass]) & is_cg
# Causal variants are drawn independently per variant, so a causal gene can come
# out with ZERO of them -- it then has h2 exactly 0 while still counting as a
# positive control, which drags the control's rank statistic toward the null and
# makes a correct pipeline look like it failed. (Before this: 7 of the top-10
# causal genes in the smoke run were h2 = 0.) Force one per empty causal gene,
# sampled with the same per-class propensity so the effect-class mix is kept.
n_forced = 0
for _gi in causal_genes:
    _m = np.where(gcode == _gi)[0]
    if causal[_m].any():
        continue
    _w = PI[eclass[_m]]
    causal[_m[rng.choice(len(_m), p=_w / _w.sum())]] = True
    n_forced += 1
_have = np.array([bool(causal[gcode == g].any()) for g in causal_genes])
assert _have.all(), f"{(~_have).sum()} causal genes still have no causal variant"
draw = rng.normal(0, SD[eclass] * scale)
# class 2 (pLoF-like) is directional: 90% positive
flip = np.where(rng.random(len(cod_idx)) < 0.90, 1.0, -1.0)
draw = np.where(eclass == 2, flip * np.abs(draw), draw)
beta_c[causal] = draw[causal]
say(f"{n_causal_genes} causal genes; {int(causal.sum()):,} causal coding variants "
    f"({causal.sum() / max(is_cg.sum(), 1) * 100:.1f}% of variants in causal genes); "
    f"{n_forced} gene(s) needed a forced variant")

def dot_chunked(mat, b, n_out, step=2000):
    """sum_j b_j * mat[j, :] without materialising mat as float64 (2.6 GB)."""
    out = np.zeros(n_out)
    for i in range(0, mat.shape[0], step):
        out += b[i:i + step] @ mat[i:i + step].astype(np.float64)
    return out


# Every causal coding variant lies in a causal gene (causal = ... & is_cg), so
# the coding genetic value is exactly the sum of the per-gene burdens. Building
# them individually costs one pass and makes the per-gene rescale below free.
burden = {int(gi): dot_chunked(D[(gcode == gi) & causal],
                               beta_c[(gcode == gi) & causal], N)
          for gi in causal_genes}
g_cod = sum(burden.values()) if burden else np.zeros(N)

# non-coding background: same spike-and-slab, no annotation, symmetric
nc = np.where(gene < 0)[0]
PI_NC = 0.001
nc_causal = nc[rng.random(len(nc)) < PI_NC]
nc_scale = (2 * maf[nc_causal] * (1 - maf[nc_causal])) ** (a.alpha / 2)
beta_nc = rng.normal(0, nc_scale)
say(f"non-coding causal variants: {len(nc_causal):,} of {len(nc):,}")


def read_bed_rows(path, rows, n_samples):
    """Read selected variant records back out of the .bed we just wrote."""
    nb = (n_samples + 3) // 4
    # plink1 code -> ALT dosage. MUST be the exact inverse of CODE above
    # (dosage 0,1,2 -> code 3,2,0), i.e. code 0=hom-ALT=2, 2=het=1, 3=hom-REF=0.
    # Getting this backwards reads 2-x and silently negates every effect size.
    unpack = np.array([2, np.nan, 1, 0])
    out = np.empty((len(rows), n_samples))
    with open(path, "rb") as fh:
        for k, r in enumerate(rows):
            fh.seek(3 + int(r) * nb)
            raw = np.frombuffer(fh.read(nb), np.uint8)
            codes = np.empty((nb, 4), np.uint8)
            for s in range(4):
                codes[:, s] = (raw >> (2 * s)) & 3
            out[k] = unpack[codes.ravel()[:n_samples]]
    return out


g_non = (read_bed_rows(f"{a.outdir}/{a.prefix}.bed", nc_causal, N).T @ beta_nc
         if len(nc_causal) else np.zeros(N))
say(f"raw genetic variance: coding {g_cod.var():.4g}, non-coding {g_non.var():.4g}")

# ------------------------------------------------- 7. scale to target h2
ncp_target = np.zeros(a.n_genes)
if a.ncp_ladder:
    # Per-gene NCP ladder. A single --h2 gives every causal gene roughly the same
    # power, so at large N the positive control saturates (every method finds
    # every gene, p ~ 1e-25) and stops discriminating between BURDEN/SKAT/SKAT-O
    # exactly when it is most needed. A log-uniform ladder instead spreads genes
    # across the whole power curve, so the control becomes a rank correlation
    # over genes rather than a hit count. NCP ~ N * h2_gene for the 1-df burden
    # test (exact only as h2_total -> 0; h2_total here is ~1%).
    lo, hi = (float(x) for x in a.ncp_ladder.split(","))
    assert 0 < lo < hi, (lo, hi)
    draws = np.exp(rng.uniform(np.log(lo), np.log(hi), len(causal_genes)))
    for gi, tgt in zip(causal_genes, draws):
        gi = int(gi)
        ncp_target[gi] = tgt
        v = burden[gi].var()
        if v <= 0:
            continue
        f = np.sqrt((tgt / N) / v)
        burden[gi] *= f
        beta_c[(gcode == gi) & causal] *= f
    g_cod = sum(burden.values()) if burden else np.zeros(N)
    say(f"NCP ladder log-uniform on [{lo}, {hi}]: coding h2 is now an OUTCOME "
        f"({g_cod.var():.4g}), --h2={a.h2} sets the non-coding part only")
    c_non = np.sqrt(a.h2 * (1 - a.frac_h2_coding) / g_non.var()) if g_non.var() > 0 else 0.0
    beta_nc *= c_non
    g_non *= c_non
else:
    c_cod = np.sqrt(a.h2 * a.frac_h2_coding / g_cod.var()) if g_cod.var() > 0 else 0.0
    c_non = np.sqrt(a.h2 * (1 - a.frac_h2_coding) / g_non.var()) if g_non.var() > 0 else 0.0
    beta_c *= c_cod
    beta_nc *= c_non
    g_cod *= c_cod
    g_non *= c_non
    for gi in causal_genes:
        burden[int(gi)] *= c_cod
g = g_cod + g_non
H2 = g.var()                      # realised chr20 h2; the ladder makes it an outcome
if not a.ncp_ladder and abs(H2 - a.h2) > 0.05 * a.h2:
    # Usually NOT "a component had zero variance": g_cod and g_non are computed
    # from the same genotypes and are correlated, so var(g_cod + g_non) is not
    # the sum of the two targets.
    print(f"WARNING: realised chr20 h2 {H2:.4g} != target {a.h2}; components are "
          f"correlated (cov={np.cov(g_cod, g_non)[0, 1]:.3g}) or one had zero variance",
          flush=True)
say(f"scaled: var(g)={H2:.4g}; beta_c range [{beta_c.min():.3f}, {beta_c.max():.3f}]")

# per-gene true h2 (variance of the gene's own burden of causal effects)
gene_h2 = np.zeros(a.n_genes)
for gi in causal_genes:
    gene_h2[int(gi)] = burden[int(gi)].var()
# List CAUSAL genes only. Ranking all n_genes and taking the top 10 pads the
# table with null genes at h2 = 0 exactly whenever there are fewer than 10
# causal genes, which reads as "causal genes with no effect" and is not.
_ord = causal_genes[np.argsort(-gene_h2[causal_genes])]
_n0 = int((gene_h2[causal_genes] == 0).sum())
_bmul = a.h2_binary_mult * BIN_NCP_FRAC
say(f"causal genes by true h2; {_n0} of {len(causal_genes)} have h2 == 0. "
    f"NCP_q = N*h2 at N={N}; NCP_b = NCP_q * {_bmul:.3f}:")
say("   gene       h2         NCP_target   NCP_quant   NCP_binary")
for gi in _ord[:12]:
    gi = int(gi)
    nq = N * gene_h2[gi]
    say(f"   GENE{gi + 1:04d}  {gene_h2[gi]:.3e}  "
        f"{ncp_target[gi]:9.1f}  {nq:10.1f}  {nq * _bmul:10.1f}")
if len(_ord) > 12:
    say(f"   ... {len(_ord) - 12} more causal genes")
if a.ncp_ladder:
    _tg = ncp_target[causal_genes]
    _rl = N * gene_h2[causal_genes]
    _ok = _tg > 0
    _err = np.abs(_rl[_ok] - _tg[_ok]) / _tg[_ok]
    say(f"ladder delivered: max relative NCP error {_err.max():.2e} "
        f"over {_ok.sum()} genes (0 = exact)")

# ------------------------------------------------------- 8. phenotypes
age = np.clip(rng.normal(50, 8, N), 20, 80)
sex01 = (sex == 1).astype(float)                 # 1 = male
z_age = (age - age.mean()) / age.std()
cov_raw = 0.7 * z_age + 1.0 * (sex01 - sex01.mean())
cov = cov_raw * np.sqrt(a.var_covar / cov_raw.var())

# background polygenic effect of the (unsimulated) rest of the genome,
# transmitted through the pedigree so it lands exactly in the span of the GRM
poly = np.zeros(N)
sd_p = np.sqrt(a.h2_poly)
poly[:a.n_singletons] = rng.normal(0, sd_p, a.n_singletons)
pv = rng.normal(0, sd_p, 2 * a.n_families)
poly[row_par] = pv
mid = pv.reshape(-1, 2).mean(1)
poly[row_chi] = np.repeat(mid, 2) + rng.normal(0, sd_p / np.sqrt(2), n_child)
say(f"var(poly)={poly.var():.3f} (target {a.h2_poly})")

def make_traits(h2_target, tag):
    s = np.sqrt(h2_target / H2) if H2 > 0 else 0.0
    ve = 1.0 - h2_target - a.h2_poly - a.var_covar
    assert ve > 0, (h2_target, ve)
    liab = s * g + poly + cov + rng.normal(0, np.sqrt(ve), N)
    yq = liab
    # The binary arm gets its OWN liability when --h2-binary-mult != 1, because
    # thresholding the quantitative one would hand it 14x less non-centrality at
    # K=0.01. Same genes and same relative effects, only the genetic scale moves.
    if a.h2_binary_mult == 1.0:
        liab_b = liab
    else:
        h2b = h2_target * a.h2_binary_mult
        # The liability has a fixed budget of 1. Scaling the genetic component up
        # 14x is feasible for the ladder arm (h2 ~ 0.01) and NOT for a strong
        # positive-control arm, so clamp -- loudly, because a silently weakened
        # positive control is exactly the failure this amendment exists to avoid.
        h2b_max = 1.0 - a.h2_poly - a.var_covar - 0.05
        if h2b > h2b_max:
            print(f"WARNING: {tag}: binary liability h2 {h2b:.3g} "
                  f"(= {h2_target:.3g} x {a.h2_binary_mult:.4g}) exceeds the budget "
                  f"1 - h2_poly - var_covar - 0.05 = {h2b_max:.3g}; clamped, so the "
                  f"effective multiplier is {h2b_max / h2_target:.2f} not "
                  f"{a.h2_binary_mult:.4g} and this arm's NCP ratio is "
                  f"{h2b_max / h2_target * BIN_NCP_FRAC:.3f}x, not "
                  f"{a.h2_binary_mult * BIN_NCP_FRAC:.3f}x", flush=True)
            h2b = h2b_max
        ve_b = 1.0 - h2b - a.h2_poly - a.var_covar
        assert ve_b > 0, (h2b, ve_b)
        liab_b = (np.sqrt(h2b / H2) * g + poly + cov
                  + rng.normal(0, np.sqrt(ve_b), N))
    cut_ = np.quantile(liab_b, 1 - a.prevalence)
    yb = (liab_b > cut_).astype(int)
    say(f"{tag}: var(y)={yq.var():.3f}  cases={yb.sum()} ({yb.mean() * 100:.2f}%)")
    # A low prevalence at small N is degenerate, and it fails LOUDLY in the
    # statistics rather than in the generator: measured at N=600 / K=0.01 (6
    # cases) the binary positive-control arm came back with lambda_GC ~ 100 on
    # markers whose true effect is ZERO, under BOTH an exact pedigree GRM and an
    # estimated one. That is the case count, not the GRM and not SAIGE.
    if 0 < yb.sum() < 100:
        print(f"WARNING: {tag}: only {yb.sum()} cases at prevalence "
              f"{a.prevalence}. Binary results at this scale are not "
              f"interpretable -- expect wildly inflated test statistics even on "
              f"null markers. Raise --prevalence or N for any binary claim.",
              flush=True)
    return yq, yb


# With a ladder, the primary arm is generated at the REALISED h2, so each gene's
# drawn NCP is delivered exactly rather than rescaled by h2_target/H2.
H2_PRIMARY = H2 if a.ncp_ladder else a.h2
yq, yb = make_traits(H2_PRIMARY, f"h2={H2_PRIMARY:.4g}")
# Kept as the saturated "must be found by anything" smoke arm; with a ladder it
# is deliberately redundant with the top of the ladder.
yq_pc, yb_pc = make_traits(a.h2_pc, f"h2={a.h2_pc} (positive control)")
# negative control: same pedigree-transmitted background and covariates, but
# ZERO chr20 heritability -- step 2 p-values on it must be uniform.
y_null = poly + cov + rng.normal(0, np.sqrt(1 - a.h2_poly - a.var_covar), N)

# ------------------------------------------------------------- 9. outputs
with open(f"{a.outdir}/pheno.txt", "w") as fh:
    fh.write("FID IID y_quant y_binary y_quant_pc y_binary_pc y_null age sex\n")
    for i in range(N):
        fh.write(f"{fid[i]} {iid[i]} {yq[i]:.6f} {yb[i]} {yq_pc[i]:.6f} {yb_pc[i]} "
                 f"{y_null[i]:.6f} {age[i]:.2f} {int(sex01[i])}\n")

beta_all = np.zeros(n_kept)
beta_all[cod_idx] = beta_c
beta_all[nc_causal] = beta_nc
anno_all = np.full(n_kept, "noncoding", dtype=object)
anno_all[cod_idx] = ANN[anno]
ecl_all = np.full(n_kept, -1, np.int8)
ecl_all[cod_idx] = eclass
gname = np.where(gene >= 0, np.char.add("GENE", np.char.zfill((gene + 1).astype(str), 4)), ".")
is_causal_gene = np.isin(gene, causal_genes)

with gzip.open(f"{a.outdir}/variants.tsv.gz", "wt") as fh:
    fh.write("ID\tCHR\tPOS\tAC\tAF\tMAF\tgene\tcausal_gene\tanno\ttrue_class\tbeta\n")
    for i in range(n_kept):
        fh.write(f"{bim_ids[i]}\t20\t{pos[i]}\t{ac[i]}\t{af[i]:.6g}\t{maf[i]:.6g}\t"
                 f"{gname[i]}\t{int(is_causal_gene[i])}\t{anno_all[i]}\t{ecl_all[i]}\t"
                 f"{beta_all[i]:.6g}\n")

# SAIGE group file, keyed on bim column 2
with open(f"{a.outdir}/groupFile.txt", "w") as fh:
    order = np.argsort(gene[cod_idx], kind="stable")
    gs = gene[cod_idx][order]
    bnd = np.searchsorted(gs, np.arange(a.n_genes), "left")
    bnd2 = np.searchsorted(gs, np.arange(a.n_genes), "right")
    ngw = 0
    for gi in range(a.n_genes):
        sel = cod_idx[order[bnd[gi]:bnd2[gi]]]
        if len(sel) < 2:
            continue
        ngw += 1
        nm = f"GENE{gi + 1:04d}"
        fh.write(nm + " var " + " ".join(bim_ids[sel]) + "\n")
        fh.write(nm + " anno " + " ".join(ANN[anno[np.searchsorted(cod_idx, sel)]]) + "\n")
say(f"group file: {ngw} genes")

with open(f"{a.outdir}/truth_genes.tsv", "w") as fh:
    # NCP_quant/NCP_binary are what the Spearman(NCP, -log10 p) control is run
    # against; NCP_target is the ladder draw, so target-vs-realised is auditable.
    fh.write("gene\tcausal\tn_var\tn_high\tn_medium\tn_low\tn_causal_var\th2\t"
             "NCP_target\tNCP_quant\tNCP_binary\n")
    _bmul = a.h2_binary_mult * BIN_NCP_FRAC
    for gi in range(a.n_genes):
        m = gcode == gi
        if not m.any():
            continue
        nq = N * gene_h2[gi]
        fh.write(f"GENE{gi + 1:04d}\t{int(gi in causal_genes)}\t{int(m.sum())}\t"
                 f"{int((anno[m] == 2).sum())}\t{int((anno[m] == 1).sum())}\t"
                 f"{int((anno[m] == 0).sum())}\t{int((causal & m).sum())}\t"
                 f"{gene_h2[gi]:.4e}\t{ncp_target[gi]:.2f}\t{nq:.2f}\t"
                 f"{nq * _bmul:.2f}\n")

json.dump({"N": N, "n_variants": int(n_kept), "n_coding": int(len(cod_idx)),
           "n_genes": a.n_genes, "n_causal_genes": int(n_causal_genes),
           "alpha": a.alpha, "h2": a.h2, "h2_pc": a.h2_pc, "h2_poly": a.h2_poly,
           "frac_h2_coding": a.frac_h2_coding, "prevalence": a.prevalence,
           "seed": a.seed, "n_families": a.n_families,
           "n_singletons": a.n_singletons, "n_forced_causal_var": int(n_forced),
           "h2_realised": float(H2), "h2_primary_arm": float(H2_PRIMARY),
           "ncp_ladder": a.ncp_ladder, "h2_binary_mult": a.h2_binary_mult,
           "binary_ncp_frac": float(BIN_NCP_FRAC),
           # the genome-wide SFS is deliberately distorted by thinning; it is
           # EXACT inside genes and at MAF >= thin_maf (amendment A5)
           "sfs_distorted_outside_genes": a.thin_intergenic < 1.0,
           "thin_intergenic": a.thin_intergenic, "thin_maf": a.thin_maf,
           "n_thinned": int(n_thinned)},
          open(f"{a.outdir}/sim_meta.json", "w"), indent=2)

# n_coding * N bytes of scratch (2.5 GB at N=50,000); nothing downstream reads it
del D
try:
    os.remove(cod_path)
except OSError:
    pass
say("done")
