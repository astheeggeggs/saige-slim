#!/usr/bin/env python3
"""Write the PINNED-SAMPLE fixture (LEDGER #174): a phenotype file whose extra
covariate separates a group of controls, and two custom-weight group files in
which the separated carriers hold the heavy weights.

WHY A FIXTURE HAD TO BE BUILT. No bundled model has a fitted probability below
0.005, so no golden could see a saddlepoint term of a sample the null model
holds certain. Separation makes them: a covariate that is 1 only for samples
with no case drives their fitted values towards 0 for as long as the fit keeps
iterating.

WHAT IT IS. `sep` is 1 for SEP_N = 60 controls -- every carrier of the five
HEAVY markers (41 samples; the script REFUSES if any of them is a case), filled
up in phenotype-file order with the first controls that carry none of them --
and 0 for everyone else. 60 of 1,000 is what step 1's glm needs to take a
useful share of the group below 1e-8: R's glm stops on the relative change of
the deviance, the separated group's share of which is ~2 n mu, so the group's
fitted values stop near 1e-8 x deviance / (1.26 n) whatever the data say. With
8 separated samples they stop at 1.7e-8..4.8e-7 and NONE is below the fit
gate's 1e-8; with 60 they land at 1.7e-9..1.2e-7 and 13 are. The group
straddling the threshold is the fixture's finding, not an accident of it.

THE TWO GROUP FILES share GENE1's markers and annotations with
group_new_snpid_withWeights.txt and differ in the weight line only:

  _pinned.txt       heavy 1, LIGHT 0.01, the rest 1e-4  (contrast 100)
  _pinnedTail.txt   heavy 1, LIGHT 1e-4, the rest 1e-6  (contrast 10,000)

LIGHT is four lof markers whose ten carriers include four cases, so the cell's
score clears --SPAcutoff on the light carriers alone. Neither cell is one-atom
(41 heavy carriers: LEDGER #172's routing does not take it). The contrast
decides what the fitted model can show. The separated carriers' share of the
score variance is ~ contrast^2 x mu, so at 1e4 a group fitted to ~1e-8 swamps
the light carriers and no saddlepoint fires; that file is for the PLANTED model
(harness/tools/plant_pinned_model.R, fitted values at 1e-12), where the
saddlepoint collapses to 1e-8..1e-13 against a reference of 0.02.

GENE2 is given weight 1 throughout: a custom-weight cell with no heavy/light
contrast, which the guard must leave where it was.

    python3 harness/tools/make_pinned_fixture.py            # write the three files
    python3 harness/tools/make_pinned_fixture.py --check    # compare, write nothing
"""
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IN = os.path.join(REPO, "extdata", "input")
PHENO_IN = os.path.join(IN, "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt")
GROUP_IN = os.path.join(IN, "group_new_snpid_withWeights.txt")
BED = os.path.join(IN, "genotype_100markers")
PHENO_OUT = os.path.join(IN, "pheno_1000samples_pinned.txt")
GROUP_OUT = os.path.join(IN, "group_new_snpid_withWeights_pinned.txt")
GROUP_TAIL_OUT = os.path.join(IN, "group_new_snpid_withWeights_pinnedTail.txt")

HEAVY = ["rs31", "rs35", "rs45", "rs46", "rs49"]
LIGHT = ["rs10", "rs32", "rs34", "rs44"]
SEP_N = 60
MAGIC = b"\x6c\x1b\x01"


def carriers(markers):
    """IIDs carrying at least one A1 allele at any of `markers` (hard calls)."""
    with open(BED + ".fam") as fh:
        iids = [line.split()[1] for line in fh]
    with open(BED + ".bim") as fh:
        ids = [line.split()[1] for line in fh]
    n = len(iids)
    bpm = (n + 3) // 4
    with open(BED + ".bed", "rb") as fh:
        body = fh.read()
    if body[:3] != MAGIC:
        sys.exit("not a variant-major PLINK 1 bed: " + BED + ".bed")
    out = set()
    for m in markers:
        j = ids.index(m)
        rec = body[3 + j * bpm: 3 + (j + 1) * bpm]
        for i in range(n):
            code = (rec[i // 4] >> (2 * (i % 4))) & 3
            if code == 1:
                sys.exit("missing call at %s, sample %s: the fixture assumes none" % (m, iids[i]))
            if code in (0, 2):          # 00 hom A1, 10 het
                out.add(iids[i])
    return out


def build():
    with open(PHENO_IN) as fh:
        lines = fh.read().splitlines()
    header = lines[0].split()
    iid_col, y_col = header.index("IID"), header.index("y_binary")
    rows = [ln.split() for ln in lines[1:]]
    in_pheno = set(r[iid_col] for r in rows)
    heavy = carriers(HEAVY) & in_pheno
    cases = [r[iid_col] for r in rows if r[iid_col] in heavy and r[y_col] != "0"]
    if cases:
        sys.exit("a carrier of a heavy marker is a case (%d): the group would not be separated" % len(cases))
    sep = set(heavy)
    for r in rows:                       # fill in file order with controls carrying none of them
        if len(sep) >= SEP_N:
            break
        if r[y_col] == "0" and r[iid_col] not in sep:
            sep.add(r[iid_col])
    if len(sep) != SEP_N:
        sys.exit("separated group has %d samples, not %d" % (len(sep), SEP_N))
    pheno = [lines[0] + " sep"]
    for ln, r in zip(lines[1:], rows):
        pheno.append(ln + (" 1" if r[iid_col] in sep else " 0"))

    with open(GROUP_IN) as fh:
        g = fh.read().splitlines()

    def group(light, rest):
        out = []
        for k in range(0, len(g), 3):
            var, anno = g[k].split(), g[k + 1]
            gene, ids = var[0], var[2:]
            if gene == "GENE1":
                w = ["1" if m in HEAVY else (light if m in LIGHT else rest) for m in ids]
            else:
                w = ["1"] * len(ids)
            out += [g[k], anno, " ".join([gene, "weight"] + w)]
        return out

    return {PHENO_OUT: pheno, GROUP_OUT: group("0.01", "1e-4"), GROUP_TAIL_OUT: group("1e-4", "1e-6")}, len(heavy)


def main():
    check = len(sys.argv) > 1 and sys.argv[1] == "--check"
    files, n_heavy = build()
    bad = 0
    for path, content in files.items():
        text = "\n".join(content) + "\n"
        if check:
            same = os.path.exists(path) and open(path).read() == text
            print(("ok   " if same else "DIFF ") + os.path.relpath(path, REPO))
            bad += 0 if same else 1
        else:
            with open(path, "w") as fh:
                fh.write(text)
            print("wrote " + os.path.relpath(path, REPO))
    print("separated group: %d samples, %d of them carriers of %s" % (SEP_N, n_heavy, ",".join(HEAVY)))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
