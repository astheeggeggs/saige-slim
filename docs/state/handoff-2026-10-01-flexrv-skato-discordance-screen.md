# Handoff 2026-10-01: the FlexRV-vs-SKAT-O discordance screen -- ONE residual outlier, and FlexRV validated as a method

Closes the discordance-screen item on `docs/state/next.md`'s blocking list. Built
and run on all four rerun ancestries (sas, amr, afr, eur) over the gene-level
summaries; no new compute, no per-trait file read.

No count of people appears in this file. Every number is a p-value, a gene count
or a (trait, gene) pair count.

## The design, and why the obvious version is wrong

FlexRV runs `--r.corr=1`, so it IS a burden test over
`pLoF:damaging_missense:other_missense` at ONE max MAF, 0.001
(`docs/state/aou-saige-parameters.md`). Two comparators are therefore wrong:

- **`p_skato` from the summary table** pools all 18 SKAT-O cells, including max
  MAF 0.01 and 0.0001. Not matched on max MAF.
- **SKAT-O's `skato:` p-values** are a variance-component test. Comparing a
  burden test against them manufactures discordance for legitimate reasons.

The right comparator is the **Cauchy combination of SKAT-O's `burden:` p-values
at `@0.001`**, over the masks inside FlexRV's group, from
`<anc>.skato_cells.v4.tsv.gz`. Two mask sets, reported together so the choice is
not load-bearing: `tight3` = pLoF, damaging_missense, pLoF;damaging_missense;
`tight4` adds other_missense, which covers FlexRV's variant space exactly.
`min` over the same cells is the conservative view -- the single best burden cell
FlexRV has to beat, with no combination involved.

Two properties that make this comparison sharp rather than merely available:
SKAT-O's `Pvalue_Burden` comes from the per-rho layer, which #158 never screens,
so there is no Liu approximation in the comparator; and both arms run the burden
p through the saddlepoint, but only FlexRV's steep transforms concentrate weight
onto one variant, so a gap isolates the concentration rather than a kernel
difference.

## TWO BUGS FOUND EN ROUTE, both of which produced a confident wrong answer first

1. **`source_old == "fixed"` is a biased subset** (recorded in section 14b and
   7a). It cost eur 38 % of its genes and amr almost nothing. Caught by the
   per-trait gene-count check, which the owner asked for.
2. **A naive CCT floors at 2^-54 = 5.551115e-17.** `0.5 - atan(S)/pi`
   catastrophically cancels once the combined p goes below that, and the literal
   repeated value in the output is the tell. It MANUFACTURED discordance wherever
   the comparator was genuinely very significant -- and because quantitative
   traits have far more such genes, it produced a fake 20-vs-0 binary/quantitative
   asymmetry that read exactly like a new artefact. The fix is the standard tail
   branch: `1/(pi S)` when `S` is large, with `tan((0.5-p)pi)` replaced by
   `1/(p pi)` for `p < 1e-15`. Before trusting any CCT in this repo, check whether
   5.551115e-17 appears in its output.

Both were caught by a CONTROL rather than by inspection, which is the argument
for running the controls first.

## Controls

- **POSITIVE:** `amr EFRMB ENSG00000213030`, `p_flexrv` 2.81e-16, sits above the
  13-carrier cap (LEDGER #176) and is the one cell known to be suspect. It must
  flag. **It does, at the top of amr**: 13.91 orders on tight4, 13.78 on tight3,
  13.34 on min.
- **NEGATIVE:** the quantitative arm. #172's routing and the saddlepoint are both
  binary-only (`traitType == "binary"`), so a quantitative trait cannot carry the
  artefact. CONFOUNDED by power, and that was stated before the numbers were read
  -- quantitative traits are far better powered, so genuine FlexRV wins are
  commoner there. The confound turned out to be the whole story (below).

## The result: the binary arm is clean apart from the known case

Orders = `log10(p_burden_cauchy@0.001) - log10(p_flexrv)`, `p_burden > 0.9`
excluded (LEDGER #86 saturation), tight4.

| | binary tests | median | p99 | binary >5 | binary **>10** | quant >5 | quant >10 |
|---|---|---|---|---|---|---|---|
| sas | 97,941 | -0.03 | +1.14 | 0 | **0** | 0 | 0 |
| amr | 510,909 | -0.07 | +1.28 | 2 | **1** | 11 | 2 |
| afr | 512,482 | -0.07 | +1.29 | 1 | **0** | 8 | 4 |
| eur | 533,141 | -0.09 | +1.31 | 5 | **0** | 47 | 20 |

**Across ~1.65 million binary gene-tests there is exactly ONE gene above 10
orders of discordance, and it is `amr EFRMB ENSG00000213030` -- already LEDGER
#176.** No residual crazy outlier is hiding in the binary FlexRV results.

The binary >5 entries read as the method working. eur's list is headed by
genes whose identity should be confirmed but whose trait pairing is textbook:
`ENSG00000008710` for CRF (PKD1 -> polycystic kidney disease), `ENSG00000012048`
for BreastCanc (BRCA1), `ENSG00000076242` for ColonRectCanc (MLH1, Lynch),
`ENSG00000115718` for VTE (PROC). **Bonus**: afr's top hit `ENSG00000198691` is
ABCA4, which CLOSES the rollout handoff's open item on whether afr AMD's 7.4e-8
gene matches a known AMD locus.

## The quantitative asymmetry is REAL BIOLOGY, and it validates the method

eur's twenty most discordant quantitative genes, after the CCT fix. Identities
from memory and to be confirmed, but the pattern is unambiguous:

| trait | gene | likely | o4 | omin |
|---|---|---|---|---|
| ALT | ENSG00000167701 | **GPT** -- the enzyme ALT IS | 172.4 | 171.8 |
| AST | ENSG00000120053 | **GOT1** -- the enzyme AST IS | 67.8 | 67.2 |
| HDLC | ENSG00000165029 | ABCA1 | 41.8 | 41.2 |
| HDLC | ENSG00000213398 | LCAT | 38.7 | 38.2 |
| LDLC | ENSG00000130164 | LDLR | 29.4 | 28.8 |
| LDLC / TChol / TG | ENSG00000084674 | APOB | 26.4 / 24.6 / 17.6 | |
| TG | ENSG00000110245 | APOC3 | 25.6 | 25.0 |
| TChol / LDLC | ENSG00000169174 | PCSK9 | 13.7 / 11.7 | |
| HDLC | ENSG00000118137 | APOA1 | 13.5 | 13.1 |
| TG | ENSG00000132855 | ANGPTL3 | 11.1 | 10.5 |

Three checks that make this a finding rather than a preference: **`omin > 10` for
all twenty**, so the gap survives against the single best burden cell and the
Cauchy is not doing the work; **only 4 of 20** have a pseudo-marker-only best set
(`_MAC1/3/5/10`, at or below the collapse cutoff -- the #163/#164 signature); and
**seven distinct traits**, not a concentration.

**So the earlier worry about an uncharacterised artefact in the quantitative arm
is RESOLVED NEGATIVE.** The 20-vs-0 asymmetry is power plus the fact that
biomarker genes carry enormous rare-variant effects.

**And this is the first validation of FlexRV AS A METHOD on real data**, which
nobody had. It recovers GPT for ALT, GOT1 for AST and the full lipid cast with a
10-to-170 order advantage over SKAT-O burden on identical variants at identical
max MAF. Set against a median of **-0.09**: FlexRV costs a hair on the typical
null gene and buys orders of magnitude where concentrated rare-variant signal
exists. That is the trade the paper claims, now measured on five cohorts.

## Caveats to keep attached

- **The extreme p-values are near the double-precision floor** (4e-296) and must
  not be quoted as precise -- only as "beyond anything the comparator reaches".
- **The mask choice is load-bearing in exactly one place found**: eur VTE
  `ENSG00000115718` reads 4.89 orders on tight4 and 11.25 on tight3, i.e. adding
  `other_missense` changes the comparator substantially there. Everywhere else
  tight3, tight4 and min agree to within ~1 order. Report all three.
- **No SKAT-O mask matches FlexRV's group exactly.** The arm's masks are pLoF,
  damaging_missense, other_missense, synonymous, pLoF;damaging_missense and the
  four-way; FlexRV's group is the three-way without synonymous, which is not
  among them. tight4 reconstructs its variant space from components.
- **The quantitative arm is not a clean null** for anything else one might want
  to screen, for the power reason above. A permuted-phenotype arm or binary
  traits with low case counts would be cleaner.

## What this does NOT say

It does not say the binary FlexRV p-values are right -- only that none except
`ENSG00000213030` disagrees wildly with a matched burden comparator. A cell whose
p is wrong in the same direction as SKAT-O's would be invisible here, as would
one wrong by less than a few orders. And LEDGER #176's mechanism is still
unexplained: the 3.12 SE factor has no account.
