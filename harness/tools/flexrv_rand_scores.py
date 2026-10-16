#!/usr/bin/env python3
"""Rewrite a FlexRV group file's score line with seeded random scores: the FlexRV null control.

FlexRV pools 112-192 weightings of one gene into p_FlexRV, and a production run has no null
arm of its own (synonymous variants carry no AlphaMissense score, so on the real file every
score transform would produce empty sets and only the MAF ladder would run). This tool keeps
every `var` and `anno` line byte for byte and replaces the `score:<NAME>` line by
`score:RAND` with a value per variant drawn uniformly on (0, 1), seeded by a hash of the
variant ID and SEED, so the same variant gets the same value in every file and every rerun.
Two runs use it (runbook section 19.7):
  the machinery's null     --annotation_in_groupTest=synonymous --flexRV_score=RAND
                           (p_FlexRV must be uniform: this is the synonymous control FlexRV lacks)
  the score's positive ctl --annotation_in_groupTest=<the production group> --flexRV_score=RAND
                           against the production AM run on the same genes: what the score buys
Usage: flexrv_rand_scores.py IN.txt.gz OUT.txt.gz [--seed N] [--score-name RAND]
Self-test: flexrv_rand_scores_selftest.sh.
"""
import argparse, gzip, hashlib, sys


def rand_unit(variant, seed):
    """Uniform on (0, 1), deterministic in (variant, seed): the first 8 bytes of a SHA-256."""
    h = hashlib.sha256(f"{seed}|{variant}".encode()).digest()
    u = int.from_bytes(h[:8], "big") / float(1 << 64)
    return min(max(u, 1e-6), 1 - 1e-6)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("infile"); ap.add_argument("outfile")
    ap.add_argument("--seed", type=int, default=20261016)
    ap.add_argument("--score-name", default="RAND")
    a = ap.parse_args()
    opener = gzip.open if a.infile.endswith(".gz") else open
    wopener = gzip.open if a.outfile.endswith(".gz") else open
    n_var = n_score = n_other = 0
    cur = None
    with opener(a.infile, "rt") as fi, wopener(a.outfile, "wt") as fo:
        for line in fi:
            f = line.rstrip("\n").split(" ")
            if len(f) >= 2 and f[1] == "var":
                cur = (f[0], f[2:]); n_var += 1
                fo.write(line)
            elif len(f) >= 2 and f[1].startswith("score:"):
                if cur is None or cur[0] != f[0]:
                    sys.exit(f"score line for {f[0]} without a preceding var line")
                if len(f) - 2 != len(cur[1]):
                    sys.exit(f"{f[0]}: score line has {len(f) - 2} values for {len(cur[1])} variants")
                vals = " ".join(f"{rand_unit(v, a.seed):.6f}" for v in cur[1])
                fo.write(f"{f[0]} score:{a.score_name} {vals}\n"); n_score += 1
            else:
                n_other += 1
                fo.write(line)
    print(f"{a.outfile}: {n_var} regions, {n_score} score lines replaced by score:{a.score_name} (seed {a.seed}), {n_other} other lines copied")
    if n_score != n_var:
        sys.exit(f"ERROR: {n_var} var lines but {n_score} score lines")


if __name__ == "__main__":
    main()
