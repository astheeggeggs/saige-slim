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
a = ap.parse_args()

rng = np.random.default_rng(a.seed)
t0 = time.time()


def say(*m):
    print(f"[{time.time() - t0:7.1f}s]", *m, flush=True)


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
cod_dose, cod_idx = [], []
ev_ptr = 0
n_seen = n_kept = n_multi = 0

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
    codes = CODE[dose].reshape(-1, 4)
    buf.append(codes[:, 0] | (codes[:, 1] << 2) | (codes[:, 2] << 4) | (codes[:, 3] << 6))
    if len(buf) >= BUFMAX:
        bed.write(np.concatenate(buf).tobytes())
        buf = []
    ip = int(p)
    ref, alt = var.alleles[0], var.alleles[1]
    bim.write(f"20\t20:{ip}:{ref}:{alt}\t0\t{ip}\t{alt}\t{ref}\n")
    k = np.searchsorted(ex_s, p, "right") - 1
    gi = int(ex_g[k]) if (k >= 0 and p < ex_e[k]) else -1
    pos_l.append(ip)
    ac_l.append(ac)
    gene_l.append(gi)
    if gi >= 0:
        cod_dose.append(dose.copy())
        cod_idx.append(n_kept)
    n_kept += 1

if buf:
    bed.write(np.concatenate(buf).tobytes())
bed.close()
bim.close()
say(f"sites seen {n_seen:,} (dropped {n_multi:,} non-biallelic); "
    f"kept {n_kept:,} polymorphic in EUR; coding {len(cod_idx):,}")

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

np.save(f"{a.outdir}/_pos.npy", pos)
np.save(f"{a.outdir}/_ac.npy", ac)
np.save(f"{a.outdir}/_gene.npy", gene)
bim_ids = np.array([l.split("\t")[1] for l in
                    open(f"{a.outdir}/{a.prefix}.bim")], dtype=object)
assert len(bim_ids) == n_kept
say("wrote bed/bim/fam")

# ------------------------------------------- 5. annotations + effect classes
cod_idx = np.array(cod_idx)
D = np.stack(cod_dose).astype(np.int8) if len(cod_idx) else np.zeros((0, N), np.int8)
del cod_dose
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
n_causal_genes = int(round(a.frac_causal_genes * a.n_genes))
causal_genes = np.sort(rng.choice(a.n_genes, n_causal_genes, replace=False))
is_cg = np.isin(gene[cod_idx], causal_genes)
scale = (2 * cmaf * (1 - cmaf)) ** (a.alpha / 2)      # slab sd multiplier

PI = np.array([0.02, 0.35, 0.55])                     # P(causal | true class)
SD = np.array([0.25, 1.00, 1.20])                     # slab sd  | true class
beta_c = np.zeros(len(cod_idx))
causal = (rng.random(len(cod_idx)) < PI[eclass]) & is_cg
draw = rng.normal(0, SD[eclass] * scale)
# class 2 (pLoF-like) is directional: 90% positive
flip = np.where(rng.random(len(cod_idx)) < 0.90, 1.0, -1.0)
draw = np.where(eclass == 2, flip * np.abs(draw), draw)
beta_c[causal] = draw[causal]
say(f"{n_causal_genes} causal genes; {int(causal.sum()):,} causal coding variants "
    f"({causal.sum() / max(is_cg.sum(), 1) * 100:.1f}% of variants in causal genes)")

def dot_chunked(mat, b, n_out, step=2000):
    """sum_j b_j * mat[j, :] without materialising mat as float64 (2.6 GB)."""
    out = np.zeros(n_out)
    for i in range(0, mat.shape[0], step):
        out += b[i:i + step] @ mat[i:i + step].astype(np.float64)
    return out


g_cod = dot_chunked(D, beta_c, N)

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
c_cod = np.sqrt(a.h2 * a.frac_h2_coding / g_cod.var()) if g_cod.var() > 0 else 0.0
c_non = np.sqrt(a.h2 * (1 - a.frac_h2_coding) / g_non.var()) if g_non.var() > 0 else 0.0
beta_c *= c_cod
beta_nc *= c_non
g_cod *= c_cod
g_non *= c_non
g = g_cod + g_non
if abs(g.var() - a.h2) > 0.05 * a.h2:
    print(f"WARNING: realised chr20 h2 {g.var():.4g} != target {a.h2} -- a component "
          f"had zero variance and could not be scaled", flush=True)
say(f"scaled: var(g)={g.var():.4g} (target h2 {a.h2}); "
    f"beta_c range [{beta_c.min():.3f}, {beta_c.max():.3f}]")

# per-gene true h2 (variance of the gene's own burden of causal effects)
gene_h2 = np.zeros(a.n_genes)
gcode = gene[cod_idx]
for gi in causal_genes:
    sel = (gcode == gi) & causal
    if sel.any():
        gene_h2[gi] = dot_chunked(D[sel], beta_c[sel], N).var()
top = np.argsort(-gene_h2)[:10]
say("top causal genes by true h2 (and expected 1-df NCP at N=%d):" % N)
for gi in top:
    say(f"   GENE{gi + 1:04d}  h2={gene_h2[gi]:.3e}  NCP={N * gene_h2[gi]:.1f}")

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
    s = np.sqrt(h2_target / a.h2)
    lin = s * g + poly + cov
    ve = 1.0 - h2_target - a.h2_poly - a.var_covar
    assert ve > 0, ve
    liab = lin + rng.normal(0, np.sqrt(ve), N)
    yq = liab
    cut_ = np.quantile(liab, 1 - a.prevalence)
    yb = (liab > cut_).astype(int)
    say(f"{tag}: var(y)={yq.var():.3f}  cases={yb.sum()} ({yb.mean() * 100:.1f}%)")
    return yq, yb


yq, yb = make_traits(a.h2, f"h2={a.h2}")
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
    fh.write("gene\tcausal\tn_var\tn_high\tn_medium\tn_low\tn_causal_var\th2\tNCP\n")
    for gi in range(a.n_genes):
        m = gcode == gi
        if not m.any():
            continue
        fh.write(f"GENE{gi + 1:04d}\t{int(gi in causal_genes)}\t{int(m.sum())}\t"
                 f"{int((anno[m] == 2).sum())}\t{int((anno[m] == 1).sum())}\t"
                 f"{int((anno[m] == 0).sum())}\t{int((causal & m).sum())}\t"
                 f"{gene_h2[gi]:.4e}\t{N * gene_h2[gi]:.2f}\n")

json.dump({"N": N, "n_variants": int(n_kept), "n_coding": int(len(cod_idx)),
           "n_genes": a.n_genes, "n_causal_genes": int(n_causal_genes),
           "alpha": a.alpha, "h2": a.h2, "h2_pc": a.h2_pc, "h2_poly": a.h2_poly,
           "frac_h2_coding": a.frac_h2_coding, "prevalence": a.prevalence,
           "seed": a.seed, "n_families": a.n_families},
          open(f"{a.outdir}/sim_meta.json", "w"), indent=2)
say("done")
