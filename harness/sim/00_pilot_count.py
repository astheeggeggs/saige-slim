#!/usr/bin/env python
"""Stage-0 pilot: count sites by category WITHOUT writing a .bed.

Measures the numbers that set the region length and the thinning probability:
kept sites/Mb after the biallelic + polymorphic-in-founders filter, the MAF/MAC
spectrum, and the CDS-resident fraction at a realistic gene density.
"""
import argparse, numpy as np, tskit

ap = argparse.ArgumentParser()
ap.add_argument("--trees", required=True)
ap.add_argument("--eur-pop", default="EUR")
ap.add_argument("--n-founders", type=int, required=True)
ap.add_argument("--genes-per-mb", type=float, default=8.4, help="real chr20 density")
ap.add_argument("--cds-mean-bp", type=float, default=1700.0)
a = ap.parse_args()

ts = tskit.load(a.trees)
pop = {p.metadata["name"]: p.id for p in ts.populations()}
node_pop = ts.tables.nodes.population
eur = [i.id for i in ts.individuals() if node_pop[i.nodes[0]] == pop[a.eur_pop]]
assert len(eur) >= a.n_founders, (len(eur), a.n_founders)
founders = np.array(eur[:a.n_founders])
nodes = np.concatenate([ts.individual(i).nodes for i in founders]).astype(np.int32)
L = ts.sequence_length
nhap = len(nodes)

# founders are the base population; the gene-drop that follows preserves
# frequencies in expectation, so the founder SFS is what matters for sizing.
ac = []
n_multi = 0
for var in ts.variants(samples=nodes):
    if len(var.alleles) != 2 or var.alleles[1] is None:
        n_multi += 1
        continue
    c = int(np.count_nonzero(var.genotypes))
    if c == 0 or c == nhap:
        continue
    ac.append(c)
ac = np.array(ac)
af = ac / nhap
maf = np.minimum(af, 1 - af)
mb = L / 1e6

print(f"region {L/1e6:.2f} Mb, {a.n_founders} founders ({nhap} haplotypes)")
print(f"raw sites in ts        : {ts.num_sites:,}  ({ts.num_sites/mb:,.0f}/Mb)")
print(f"dropped multiallelic   : {n_multi:,}")
print(f"KEPT (biallelic, poly) : {len(ac):,}  ({len(ac)/mb:,.0f}/Mb)")
print()
# MAF strata that matter: the group-test window and the GRM pool
bins = [(0, 1e-4, "MAF<=1e-4 (== MAC<=10, the collapse set)"),
        (1e-4, 1e-3, "1e-4 - 1e-3"),
        (1e-3, 1e-2, "1e-3 - 1e-2"),
        (1e-2, 0.05, "1e-2 - 0.05"),
        (0.05, 0.5, "MAF>=0.05 (clean GRM pool)")]
for lo, hi, lab in bins:
    n = int(((maf > lo) & (maf <= hi)).sum())
    print(f"  {lab:<42} {n:>9,}  {n/len(ac)*100:5.2f}%   {n/mb:>8,.0f}/Mb")
n_common = int((maf >= 0.01).sum())
print(f"\nMAF>=0.01 (createSparseGRM default floor): {n_common:,} = {n_common/mb:,.0f}/Mb")
print(f"MAF>=0.05                                : {int((maf>=0.05).sum()):,} = {int((maf>=0.05).sum())/mb:,.0f}/Mb")
# categorical variance-ratio bins (default cateVarRatioMinMACVecExclude 10,20.5)
for lo, hi, lab in [(10, 20, "MAC 11-20"), (20, 10**9, "MAC>20")]:
    n = int(((ac > lo) & (ac <= hi) | ((nhap - ac > lo) & (nhap - ac <= hi))).sum()) if hi < 10**8 else int((np.minimum(ac, nhap-ac) > lo).sum())
    print(f"VR bin {lab:<10}: {n:,} ({n/mb:,.0f}/Mb)")
print()
cds_frac = a.genes_per_mb * a.cds_mean_bp / 1e6
print(f"gene model: {a.genes_per_mb}/Mb x {a.cds_mean_bp:.0f} bp CDS = {cds_frac*100:.2f}% of sequence")
print(f"  -> coding sites ~ {len(ac)/mb*cds_frac:,.0f}/Mb, {len(ac)/mb*cds_frac/a.genes_per_mb:,.0f} per gene")
