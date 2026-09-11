#!/usr/bin/env python3
"""Build the --anno-table (class, maf_lo, maf_hi, count) for 02_build_dataset.py
from real exome summary statistics. Two modes:

  --group <SAIGE region results> --n <sample size>
      Uses the per-(gene, mask, max_MAF) rows' Number_rare + Number_ultra_rare.
      Nested max_MAF cells are differenced into MAF bins, and the ultra-rare
      count splits the lowest bin at MAC <= --collapse (MAF <= collapse / 2N).
      Masks are mapped to classes with --map (default
      "pLoF=high,missense=medium,synonymous=low"; unmapped masks are ignored;
      a compound mask like "pLoF;missense" is ignored unless mapped).
  --variants <TSV> --class-col <col> --ac-col <col> --n <N> [--af-col <col>]
      Per-variant annotation joined to allele counts (or AF): direct binning
      on --bins (default "0,1e-4,1e-3,1e-2,0.5").
  --groupfile <SAIGE group file> --single <.singleAssoc.txt> [--n N]
      The SAIGE group file supplies every variant's class (its anno label,
      mapped with --map); the single-in-group file supplies AF for the
      variants that survived the ultra-rare collapse (MAC > --collapse).
      Present variants are binned exactly on --bins above collapse/2N;
      ABSENT variants are the ultra-rare set (plus anything monomorphic or
      QC-dropped -- ASSUMED to share its class mix, printed as a trend check)
      and form the [0, collapse/2N) bin. N defaults to the file's N column.

Counts are pooled over genes: the builder uses P(class | MAF bin), so the
table's N need not equal the fixture's. Prints the table and the singleton
fraction where it can be computed (variants mode only). Provenance (file
names, md5, N) goes into the table header as comment lines.
"""
import argparse, csv, hashlib, sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--group")
ap.add_argument("--variants")
ap.add_argument("--class-col", default="anno")
ap.add_argument("--ac-col", default="AC")
ap.add_argument("--af-col", default="")
ap.add_argument("--groupfile")
ap.add_argument("--anno", help="per-variant annotation table (alternative to --groupfile)")
ap.add_argument("--id-col", default="ID")
ap.add_argument("--single")
ap.add_argument("--n", type=float, default=0.0, help="sample size of the summary statistics")
ap.add_argument("--collapse", type=float, default=10.0)
ap.add_argument("--map", default="pLoF=high,missense=medium,synonymous=low")
ap.add_argument("--bins", default="0,1e-4,1e-3,1e-2,0.5")
ap.add_argument("--out", required=True)
a = ap.parse_args()
cmap = dict(kv.split("=") for kv in a.map.split(","))
classes = ("low", "medium", "high")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


rows = []                      # (class, lo, hi, count)
prov = []
if a.group:
    prov.append(f"# source: SAIGE group results {a.group} md5={md5(a.group)} N={a.n:g} collapse={a.collapse:g}")
    per = {}                   # (gene, class) -> {maxmaf: (rare, ur)}
    with open(a.group) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            if r["Group"] == "Cauchy" or r["Group"] not in cmap:
                continue
            key = (r["Region"], cmap[r["Group"]])
            per.setdefault(key, {})[float(r["max_MAF"])] = (
                int(float(r["Number_rare"])), int(float(r["Number_ultra_rare"])))
    mafs = sorted({m for d in per.values() for m in d})
    ur_hi = a.collapse / (2 * a.n)
    counts = {}
    for (gene, cls), d in per.items():
        prev_total = 0
        for k, m in enumerate(mafs):
            if m not in d:
                continue
            rare, ur = d[m]
            total = rare + ur
            if k == 0:
                # the lowest cell: ultra-rare (MAC <= collapse) and the rest
                counts[(cls, 0.0, ur_hi)] = counts.get((cls, 0.0, ur_hi), 0) + ur
                counts[(cls, ur_hi, m)] = counts.get((cls, ur_hi, m), 0) + rare
            else:
                lo = mafs[k - 1]
                counts[(cls, lo, m)] = counts.get((cls, lo, m), 0) + max(total - prev_total, 0)
            prev_total = total
    for (cls, lo, hi), c in counts.items():
        rows.append((cls, lo, hi, c))
elif a.variants:
    prov.append(f"# source: variants {a.variants} md5={md5(a.variants)} N={a.n:g} class-col={a.class_col} ac-col={a.ac_col}")
    edges = [float(x) for x in a.bins.split(",")]
    counts = {}
    n_single = {c: 0 for c in classes}
    n_all = {c: 0 for c in classes}
    with open(a.variants) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            cls = cmap.get(r[a.class_col], r[a.class_col] if r[a.class_col] in classes else None)
            if cls is None:
                continue
            if a.af_col:
                af = float(r[a.af_col])
            else:
                af = float(r[a.ac_col]) / (2 * a.n)
            maf = min(af, 1 - af)
            if maf <= 0:
                continue
            k = int(np.searchsorted(edges, maf, side="right")) - 1
            if k < 0 or k >= len(edges) - 1:
                continue
            counts[(cls, edges[k], edges[k + 1])] = counts.get((cls, edges[k], edges[k + 1]), 0) + 1
            n_all[cls] += 1
            if not a.af_col and float(r[a.ac_col]) == 1:
                n_single[cls] += 1
    for (cls, lo, hi), c in counts.items():
        rows.append((cls, lo, hi, c))
    for c in classes:
        if n_all[c]:
            print(f"singleton fraction {c}: {n_single[c] / n_all[c]:.3f} (n={n_all[c]})")
elif (a.groupfile or a.anno) and a.single:
    import gzip
    def opn(path):
        return gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    cls_of = {}
    lines = []
    if a.anno:
        with opn(a.anno) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                lab = r[a.class_col]
                if lab in cmap:
                    cls_of.setdefault(r[a.id_col], cmap[lab])
    else:
        with open(a.groupfile) as fh:
            lines = fh.read().splitlines()
    i = 0
    while i < len(lines):
        v = lines[i].split(); an = lines[i + 1].split()
        assert v[1] == "var" and an[1] == "anno" and v[0] == an[0], (i, v[:2], an[:2])
        for vid, lab in zip(v[2:], an[2:]):
            cls_of.setdefault(vid, cmap.get(lab))       # first mask wins for a repeated ID
        i += 2
        while i < len(lines) and lines[i].split()[1] not in ("var",):
            i += 1                                     # skip weight/score lines
    maf_of = {}
    n_col = []
    with opn(a.single) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            if r["MarkerID"] not in cls_of:
                continue
            af = float(r["AF_Allele2"])
            maf_of[r["MarkerID"]] = min(af, 1 - af)
            if len(n_col) < 1000:
                if "N" in r:
                    n_col.append(float(r["N"]))
                elif "N_case" in r and "N_ctrl" in r:     # a binary trait's file
                    n_col.append(float(r["N_case"]) + float(r["N_ctrl"]))
    if not a.n:
        assert n_col, "give --n: the single-variant file has no N column"
        a.n = float(np.median(n_col))
    prov.append(f"# source: annotation {a.anno or a.groupfile} md5={md5(a.anno or a.groupfile)}; single-in-group "
                f"{a.single} md5={md5(a.single)}; N={a.n:g} collapse={a.collapse:g}")
    ur_hi = a.collapse / (2 * a.n)
    edges = [e for e in [float(x) for x in a.bins.split(",")] if e > ur_hi]
    edges = [ur_hi] + edges
    counts = {}
    n_present = n_absent = n_unmapped = 0
    for vid, cls in cls_of.items():
        if cls is None:
            n_unmapped += 1
            continue
        if vid in maf_of:
            n_present += 1
            maf = maf_of[vid]
            k = int(np.searchsorted(edges, maf, side="right")) - 1
            if k < 0:
                lo, hi = 0.0, ur_hi          # present but at MAC <= collapse (a smaller N here)
            elif k >= len(edges) - 1:
                continue
            else:
                lo, hi = edges[k], edges[k + 1]
        else:
            n_absent += 1
            lo, hi = 0.0, ur_hi
        counts[(cls, lo, hi)] = counts.get((cls, lo, hi), 0) + 1
    for (cls, lo, hi), c in counts.items():
        rows.append((cls, lo, hi, c))
    n_common = sum(1 for vid, m in maf_of.items() if m > edges[-1])
    print(f"present variants above the largest bin edge {edges[-1]:g} (dropped): {n_common}")
    print(f"group-file variants: {len(cls_of)} (unmapped labels skipped: {n_unmapped}); "
          f"present in single-in-group: {n_present}; absent (ultra-rare bin, ASSUMED): {n_absent}")
    print("trend check: the high (pLoF) fraction should rise as MAF falls, and the absent bin "
          "should continue the trend of the exact bins above it")
else:
    sys.exit("give --group, --variants, or --groupfile with --single")

# drop degenerate bins (lo == hi when collapse / 2N coincides with a max_MAF
# cell) and empty ones: the builder needs a positive total per bin
rows = [r for r in rows if r[2] > r[1] and r[3] > 0]
rows.sort(key=lambda t: (t[1], t[2], classes.index(t[0])))
with open(a.out, "w") as fh:
    for p in prov:
        fh.write(p + "\n")
    fh.write("class\tmaf_lo\tmaf_hi\tcount\n")
    for cls, lo, hi, c in rows:
        fh.write(f"{cls}\t{lo:g}\t{hi:g}\t{c}\n")
tot = {}
for cls, lo, hi, c in rows:
    tot[(lo, hi)] = tot.get((lo, hi), 0) + c
print("class x MAF bin, proportions within bin (n):")
for (lo, hi) in sorted(tot):
    parts = []
    for c in classes:
        n = next((r[3] for r in rows if r[0] == c and r[1] == lo and r[2] == hi), 0)
        parts.append(f"{c}={n / tot[(lo, hi)]:.3f}" if tot[(lo, hi)] else f"{c}=nan")
    print(f"  [{lo:g}, {hi:g})  n={tot[(lo, hi)]:8d}  " + "  ".join(parts))
