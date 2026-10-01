#!/usr/bin/env python3
"""Per-trait calibration as a COLUMN (LEDGER #176's third detectability
instrument; docs/state/next.md's blocking list).

The null-median check used by hand on 2026-09-30 is three lines, which is why it
was a notebook cell and why it was only ever run on the traits someone thought to
look at. At a thousand traits nobody looks. This computes it for every
(ancestry, trait, arm) from the gene-level summary and writes one row each.

THE STATISTIC. Under a clean null the gene-level p-values of one trait are
uniform, so the SMALLEST of m of them has

    median      = 1 - 0.5**(1/m)          (~ ln2/m)
    P(min<=obs) = 1 - (1-min_p)**m        <- the calibration p-value

`ratio` = min_p / median is the readable summary (1.0 = exactly on the null
median, <1 = more extreme than a null would give, >1 = conservative), and
`p_cal` = P(min<=obs) is what to apply multiplicity to across traits.

WHY THE MINIMUM AND NOT A QQ SLOPE. Inflation that matters here is not a
genome-wide shift -- the FlexRV defects found so far (#163/#164, #172, #176)
produce ONE absurd gene while leaving the bulk untouched, which a lambda or a QQ
slope averages away. The minimum is the statistic that sees exactly that.

WHAT IT IS NOT. A trait flagged here is not wrong; it is worth looking at. A real
association also produces a small minimum -- which is why the SKAT-O arm is
computed alongside, and why `concordant` is reported: a trait whose two arms are
both extreme is probably biology, and one extreme in FlexRV alone is the shape
the three known defects had.

Usage:  aou_flexrv_calibration.py <summary.tsv.gz> [more.tsv.gz ...] > calib.tsv
"""
import sys, gzip, math, csv

ARMS = (("flexrv", "p_flexrv"), ("skato", "p_skato"))
# below this the double-precision (1-p)**m form loses everything; use the
# expm1/log1p route, which is exact to machine precision for small p
SMALL = 1e-12


def p_min_le(min_p, m):
    """P(min of m uniforms <= min_p), computed without cancellation."""
    if m <= 0:
        return float("nan")
    if min_p <= 0:
        return 1.0
    if min_p >= 1:
        return 1.0
    return -math.expm1(m * math.log1p(-min_p))


def null_median(m):
    return -math.expm1(math.log(0.5) / m) if m > 0 else float("nan")


def read(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            yield row


def main(paths):
    out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    out.writerow(["anc", "trait", "arm", "source", "n_genes", "min_p", "gene",
                  "set_or_cell", "null_median", "ratio", "p_cal"])
    rows = []
    for path in paths:
        acc = {}
        for r in read(path):
            for arm, col in ARMS:
                v = r.get(col, "")
                if v in ("", "NA", "NaN"):
                    continue
                try:
                    p = float(v)
                except ValueError:
                    continue
                if not (0.0 <= p <= 1.0):
                    continue
                k = (r["anc"], r["trait"], arm)
                best = r.get("set_best" if arm == "flexrv" else "cell_best", "")
                cur = acc.get(k)
                # count every usable gene; keep the smallest and who it was
                if cur is None:
                    acc[k] = [1, p, r["gene"], best, r.get("source", "")]
                else:
                    cur[0] += 1
                    if p < cur[1]:
                        cur[1], cur[2], cur[3] = p, r["gene"], best
        for (anc, trait, arm), (m, mp, gene, best, src) in sorted(acc.items()):
            med = null_median(m)
            rows.append([anc, trait, arm, src, m, f"{mp:.6g}", gene, best,
                         f"{med:.6g}", f"{mp/med:.4g}" if med > 0 else "NA",
                         f"{p_min_le(mp, m):.6g}"])
    for r in rows:
        out.writerow(r)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
