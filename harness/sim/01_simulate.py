#!/usr/bin/env python
"""Coalescent simulation of chr20 under a published human demography.

Species/contig/map/model all come from stdpopsim so the parameters are the
published ones rather than hand-typed:
  - HomSap chr20 (GRCh38, 64.44 Mb)
  - HapMapII_GRCh38 recombination map for chr20 (real hotspot structure)
  - OutOfAfrica_2T12 (Tennessen et al. 2012): AFR/EUR with the recent
    explosive growth (present-day EUR Ne ~5e5). Gutenkunst OutOfAfrica_3G09
    was tried first and rejected: its present-day CEU Ne is 29,725, so a
    9,000-diploid sample is ~30% of the effective population and the cohort is
    riddled with cryptic relatedness (3.2% of "unrelated" pairs exceed GRM
    0.05, skew +3.2), which makes a 0.05-cutoff sparse GRM 2.6% dense. The
    Tennessen growth also gives the realistic excess of rare variants that a
    rare-variant burden simulation depends on.

The mutation rate is taken from the DEMOGRAPHIC MODEL, not from the contig
default (1.29e-8): the Ne values were inferred under the model's own rate, so
pairing them with the lower contig rate would deflate diversity ~1.8x.

EUR is oversampled because the downstream step makes 500 nuclear families by
gene-drop, which needs 2 founders per family drawn from this pool.
"""
import argparse
import time

import stdpopsim

ap = argparse.ArgumentParser()
ap.add_argument("--left", type=float, default=None, help="window start bp")
ap.add_argument("--right", type=float, default=None, help="window end bp (truncate chr20 here)")
ap.add_argument("--n-eur", type=int, default=9000)
ap.add_argument("--model", default="OutOfAfrica_2T12")
ap.add_argument("--eur-pop", default="EUR")
ap.add_argument("--other", default="AFR:1000",
                help="comma-separated POP:n for the non-EUR populations")
ap.add_argument("--seed", type=int, default=2026)
ap.add_argument("--out", required=True)
a = ap.parse_args()

species = stdpopsim.get_species("HomSap")
model = species.get_demographic_model(a.model)
kw = {"genetic_map": "HapMapII_GRCh38", "mutation_rate": model.mutation_rate}
if a.left is not None:
    kw["left"] = a.left
if a.right is not None:
    kw["right"] = a.right
contig = species.get_contig("chr20", **kw)
samples = {a.eur_pop: a.n_eur}
for tok in filter(None, a.other.split(",")):
    k, v = tok.split(":")
    samples[k] = int(v)

print(f"contig [{contig.origin}] len={contig.length:,.0f} mu={contig.mutation_rate:.3g} "
      f"mean_r={contig.recombination_map.mean_rate:.3g}", flush=True)
print(f"model={model.id} samples={samples} seed={a.seed}", flush=True)

t0 = time.time()
ts = stdpopsim.get_engine("msprime").simulate(model, contig, samples, seed=a.seed)
print(f"simulate: {time.time() - t0:.1f}s  trees={ts.num_trees:,} "
      f"sites={ts.num_sites:,} edges={ts.num_edges:,} inds={ts.num_individuals:,}", flush=True)
ts.dump(a.out)
print(f"wrote {a.out}", flush=True)
