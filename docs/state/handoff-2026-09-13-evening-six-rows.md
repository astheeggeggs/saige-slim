# Six ledger rows fixed, three left, and where the three stand (2026-09-13, evening)

One commit per row, each gated at `both` modes locally and again by the
pre-push hook, CI green on `slim` after each push (runs 34770616125 and
34772933326). The ledger went from 9 open to **3**. Read the rows
(`python3 harness/tools/ledger.py <n>`); each STATUS cell is the authoritative
story. This file is the map.

| commit | row | one line |
|---|---|---|
| `53d77dd2` | #141 | a repeated IID in the fam/psam is refused at every step, both readers |
| `b4437000` | #140 | S9 carried to the region path: byte counts for BOTH files on the index line, flush before measuring, truncate on restart |
| `f5be6d59` | #116 | diagnosed (ingest QC gate, not the collapse) then DECIDED: header at start and one `Region Cauchy NA` row per untested region, both paths |
| `49274ad1` | #137 | the manifest refuses `--varianceRatioFile` as its help promised |
| `742726b3` | #138 | a FlexRV run refuses any set composed beside its transform sets |
| `51a8e7f2` | #136 | six help sentences corrected in every copy (a seventh twin found), six `docdrift.py` rules |

None moved a golden, each as predicted before its gate.

## What a later session must not re-learn

- **Red drills ran against a scratch install of unpatched HEAD** built by
  `git archive HEAD` + a symlinked `.pixi` + the tree's `src/*.o` copied in
  (minus the changed file's) + `R CMD INSTALL -l <scratch>`. About a minute.
  Every drill printed `system.file(package = "SAIGE")` and it was checked;
  for a C++ guard the new message string was shown ABSENT from the scratch
  `.so` and PRESENT in the patched one (`strings | grep`) before either
  colour was believed. Recipe and trap in `next.md`'s 2026-09-13 evening entry.
- **`fitNULLGLMM` pre-creates the model `.rda` EMPTY** before reading anything
  (`file.create`, `R/SAIGE_fitGLMM_fast.R:851`): a refusal test asserts SIZE 0,
  not absence.
- **The comment linter parses a comment whose first token is a bare symbol as
  code** (`# LEDGER #140 FIXED. ...` is `LEDGER` + a comment). Start such
  lines with a phrase. The same linter flags `--` inside a comment.
- **Three characterizations of #116's old behaviour named no row** and were
  found only by the gate (`test-step2-file-integrity.R`, `test-step2-region-chunkio.R`,
  `test-unit-region-degenerate.R`). Before inverting a behaviour, grep the
  tests for the SYMPTOM (`readLines(out) == 0`, `expect_null(r$pval.Region)`),
  not just the row number.
- **`R/RegionSchema.R` now states the SKAT-O results schema once** and
  `regionAssertSchema` stops any run whose rows disagree with it, so every
  committed SKAT-O golden gates the header. Add a column there first.
- **The region index line is now `... chunk N single_bytes S bytes B`**; the
  writer `indexChunkLine()` sits beside the reader `parseIndexChunkLine()` in
  `R/Util.R`. One older test pinned the legacy line verbatim and was loosened.
- **Every fastTest region run leaves `<out>.singleAssoc.txt_temp` behind.** Not
  a row, predates this session, noted here so it is not mistaken for #116.
- **`man/` is regenerated with `roxygen2::roxygenise(".", roclets = "rd")`**
  from the pinned env; it also re-wraps one unrelated paragraph (whitespace).
- **The docdrift drill at `4c2ca05a^` now reports 42 hits** (10 before the six
  #136 rules); the number is in CLAUDE.md and must move with the rule count.

## The three open rows, investigated but not fixed

**#142 -- FlexRV MAC rungs at or above `2*N*maxMAF` duplicate `unw` and enter
the pooled row.** A DECISION, not a fix: the absolute ladder is a recorded
owner's decision (`docs/design/flexrv-plan.md:236-241`). Four pooling rules are
laid out and measured by `harness/sim/flexrv_study/pool_rules.R` (new, this
session) over the study's 192-set outputs, with no SAIGE re-run:
`all` (shipped: equal-weight CCT over 192), `dedup` (drop value-identical
sets within a score transform, the row's candidate b), `apriori` (drop the
rungs at or above the bound before pooling -- computable from flags + N, so it
could be done at spec time and would also cut 25% of the work), `axis`
(weighted CCT: each score transform gets 1/16, shared over its DISTINCT
MAF-axis sets, so the score axis is not re-weighted by the duplicate count).
It reports null-arm ECDFs, power-arm detections and rank agreement per rule,
plus a control that `all` reproduces the shipped Cauchy row. Its first run's
output is appended below. What the decision turns on: at the paper's N
(~400k) no rung is degenerate, so `all` IS the paper's statistic; below
N = 250k it is not, and the sign of the distortion is not fixed (GENE1 1.45x,
GENE2 0.966x on the bundled fixture). The cheapest non-arithmetic fix is the
row's candidate (a), an N-aware warning naming the degenerate rungs, which
moves nothing; `apriori` is the cleanest arithmetic fix and moves the FlexRV
goldens.

**#139 -- chrX male doubling counts the missing-call sentinel.** Needs a
fixture: `genotype_100markers_missing` (60 of 100 markers carry missing
calls) with the bim's chromosome rewritten to X and a male sample file, the
recipe `test-step2-xnonpar-chrom.R` already uses for chrX minus the missing
fixture. Two readers, two wrong answers (PLINK: rows vanish; PGEN: a
`BETA=NA SE=NA` row beside `nan` AC/AF), plus a third path under `--condition`
that conditions on a garbage vector. Fix inside `applyMaleXnonPAR` (count after
imputation or over non-missing entries only), correct the false `cat()` at
`R/SAIGE_Test_main.R:544`, and gate with a new golden through BOTH readers
whose point is that they agree. Not started.

**#135 -- the BURDEN emitter folds a cell's AGGREGATE count as a frequency.**
Inert on every golden; bites when a cell's aggregate count approaches 2N. The
bundled 100 markers at max MAF 0.5 in ONE cell can reach it (P65's
configuration pushed further); the characterization's non-vacuity assertion is
that a cell with MORE alleles is not gated as though it had FEWER
(`--minGroupMAC_in_BurdenTest` raised above the folded value must NOT refuse
it). Fix is `MAC = altCounts` in the emitter, which is what the printed `MAC`
already is. Not started; the aggregate-count measurement was queued when the
handoff was requested.

## pool_rules.R first run (appended verbatim when it finished)

```
CONTROL: recomputed equal-weight CCT vs the file's Cauchy row -- max |dlog10 p| = 0.000475 
duplicate sets per region: value-identical median 69  a-priori rungs >= bound: 48 45 

NULL ARMS -- ECDF at alpha per rule (binomial 95% CI)
  y_null_bin_informative_flex192_c0          alpha=0.05  all=0.0429  dedup=0.0524  apriori=0.0429  axis=0.0429   [CI 0.0231, 0.0858] n=210
  y_null_bin_informative_flex192_c0          alpha=0.01  all=0.0143  dedup=0.0190  apriori=0.0143  axis=0.0190   [CI 0.0012, 0.0340] n=210
  y_null_bin_informative_flex192_c10         alpha=0.05  all=0.0429  dedup=0.0524  apriori=0.0429  axis=0.0429   [CI 0.0231, 0.0858] n=210
  y_null_bin_informative_flex192_c10         alpha=0.01  all=0.0143  dedup=0.0190  apriori=0.0143  axis=0.0190   [CI 0.0012, 0.0340] n=210
  y_null_bin_realistic_flex192_c0            alpha=0.05  all=0.0429  dedup=0.0429  apriori=0.0429  axis=0.0429   [CI 0.0231, 0.0858] n=210
  y_null_bin_realistic_flex192_c0            alpha=0.01  all=0.0095  dedup=0.0143  apriori=0.0095  axis=0.0143   [CI 0.0012, 0.0340] n=210
  y_null_bin_realistic_flex192_c10           alpha=0.05  all=0.0429  dedup=0.0429  apriori=0.0429  axis=0.0429   [CI 0.0231, 0.0858] n=210
  y_null_bin_realistic_flex192_c10           alpha=0.01  all=0.0095  dedup=0.0143  apriori=0.0095  axis=0.0143   [CI 0.0012, 0.0340] n=210
  y_null_bin_synonymous_flex192_c0           alpha=0.05  all=0.0619  dedup=0.0619  apriori=0.0619  axis=0.0619   [CI 0.0231, 0.0858] n=210
  y_null_bin_synonymous_flex192_c0           alpha=0.01  all=0.0429  dedup=0.0381  apriori=0.0381  axis=0.0381   [CI 0.0012, 0.0340] n=210
  y_null_bin_synonymous_flex192_c10          alpha=0.05  all=0.0619  dedup=0.0619  apriori=0.0619  axis=0.0619   [CI 0.0231, 0.0858] n=210
  y_null_bin_synonymous_flex192_c10          alpha=0.01  all=0.0429  dedup=0.0381  apriori=0.0381  axis=0.0381   [CI 0.0012, 0.0340] n=210
  y_null_bin_uninformative_flex192_c0        alpha=0.05  all=0.0381  dedup=0.0381  apriori=0.0381  axis=0.0381   [CI 0.0231, 0.0858] n=210
  y_null_bin_uninformative_flex192_c0        alpha=0.01  all=0.0143  dedup=0.0190  apriori=0.0143  axis=0.0190   [CI 0.0012, 0.0340] n=210
  y_null_bin_uninformative_flex192_c10       alpha=0.05  all=0.0381  dedup=0.0381  apriori=0.0381  axis=0.0381   [CI 0.0231, 0.0858] n=210
  y_null_bin_uninformative_flex192_c10       alpha=0.01  all=0.0143  dedup=0.0190  apriori=0.0143  axis=0.0190   [CI 0.0012, 0.0340] n=210
  y_null_informative_flex192_c0              alpha=0.05  all=0.0810  dedup=0.0952  apriori=0.0952  axis=0.0952   [CI 0.0231, 0.0858] n=210
  y_null_informative_flex192_c0              alpha=0.01  all=0.0190  dedup=0.0143  apriori=0.0143  axis=0.0143   [CI 0.0012, 0.0340] n=210
  y_null_informative_flex192_c10             alpha=0.05  all=0.0810  dedup=0.0952  apriori=0.0952  axis=0.0952   [CI 0.0231, 0.0858] n=210
  y_null_informative_flex192_c10             alpha=0.01  all=0.0190  dedup=0.0143  apriori=0.0143  axis=0.0143   [CI 0.0012, 0.0340] n=210
  y_null_realistic_flex192_c0                alpha=0.05  all=0.0810  dedup=0.1000  apriori=0.0905  axis=0.0952   [CI 0.0231, 0.0858] n=210
  y_null_realistic_flex192_c0                alpha=0.01  all=0.0143  dedup=0.0143  apriori=0.0143  axis=0.0143   [CI 0.0012, 0.0340] n=210
  y_null_realistic_flex192_c10               alpha=0.05  all=0.0810  dedup=0.1000  apriori=0.0905  axis=0.0952   [CI 0.0231, 0.0858] n=210
  y_null_realistic_flex192_c10               alpha=0.01  all=0.0143  dedup=0.0143  apriori=0.0143  axis=0.0143   [CI 0.0012, 0.0340] n=210
  y_null_synonymous_flex192_c0               alpha=0.05  all=0.0524  dedup=0.0524  apriori=0.0524  axis=0.0524   [CI 0.0231, 0.0858] n=210
  y_null_synonymous_flex192_c0               alpha=0.01  all=0.0000  dedup=0.0048  apriori=0.0048  axis=0.0048   [CI 0.0012, 0.0340] n=210
  y_null_synonymous_flex192_c10              alpha=0.05  all=0.0524  dedup=0.0524  apriori=0.0524  axis=0.0524   [CI 0.0231, 0.0858] n=210
  y_null_synonymous_flex192_c10              alpha=0.01  all=0.0000  dedup=0.0048  apriori=0.0048  axis=0.0048   [CI 0.0012, 0.0340] n=210
  y_null_uninformative_flex192_c0            alpha=0.05  all=0.0905  dedup=0.0905  apriori=0.0905  axis=0.0810   [CI 0.0231, 0.0858] n=210
  y_null_uninformative_flex192_c0            alpha=0.01  all=0.0238  dedup=0.0143  apriori=0.0143  axis=0.0143   [CI 0.0012, 0.0340] n=210
  y_null_uninformative_flex192_c10           alpha=0.05  all=0.0905  dedup=0.0905  apriori=0.0905  axis=0.0810   [CI 0.0231, 0.0858] n=210
  y_null_uninformative_flex192_c10           alpha=0.01  all=0.0238  dedup=0.0143  apriori=0.0143  axis=0.0143   [CI 0.0012, 0.0340] n=210

POWER ARMS -- detections at Bonferroni 0.05/n among causal genes, and Spearman rho vs 'all'
  y_binary_informative_flex192_c0            causal  63/210  all=11  dedup=11  apriori=11  axis=11   rho(dedup,apriori,axis)=0.987,0.990,0.986
  y_binary_informative_flex192_c10           causal  63/210  all=11  dedup=11  apriori=11  axis=11   rho(dedup,apriori,axis)=0.988,0.991,0.987
  y_binary_realistic_flex192_c0              causal  63/210  all=11  dedup=10  apriori=11  axis=11   rho(dedup,apriori,axis)=0.964,0.976,0.963
  y_binary_realistic_flex192_c10             causal  63/210  all=11  dedup=10  apriori=11  axis=11   rho(dedup,apriori,axis)=0.964,0.976,0.963
  y_binary_synonymous_flex192_c0             causal  63/210  all= 7  dedup= 6  apriori= 6  axis= 6   rho(dedup,apriori,axis)=0.972,0.978,0.972
  y_binary_synonymous_flex192_c10            causal  63/210  all= 7  dedup= 6  apriori= 6  axis= 6   rho(dedup,apriori,axis)=0.973,0.978,0.973
  y_binary_uninformative_flex192_c0          causal  63/210  all=10  dedup=10  apriori=10  axis=10   rho(dedup,apriori,axis)=0.973,0.988,0.974
  y_binary_uninformative_flex192_c10         causal  63/210  all=10  dedup=10  apriori=10  axis=10   rho(dedup,apriori,axis)=0.969,0.984,0.969
  y_quant_informative_flex192_c0             causal  63/210  all= 2  dedup= 2  apriori= 2  axis= 2   rho(dedup,apriori,axis)=0.966,0.982,0.965
  y_quant_informative_flex192_c10            causal  63/210  all= 2  dedup= 2  apriori= 2  axis= 2   rho(dedup,apriori,axis)=0.966,0.982,0.965
  y_quant_realistic_flex192_c0               causal  63/210  all= 2  dedup= 2  apriori= 2  axis= 2   rho(dedup,apriori,axis)=0.958,0.978,0.954
  y_quant_realistic_flex192_c10              causal  63/210  all= 2  dedup= 2  apriori= 2  axis= 2   rho(dedup,apriori,axis)=0.958,0.978,0.953
  y_quant_synonymous_flex192_c0              causal  63/210  all= 0  dedup= 0  apriori= 0  axis= 0   rho(dedup,apriori,axis)=0.970,0.977,0.970
  y_quant_synonymous_flex192_c10             causal  63/210  all= 0  dedup= 0  apriori= 0  axis= 0   rho(dedup,apriori,axis)=0.975,0.980,0.975
  y_quant_uninformative_flex192_c0           causal  63/210  all= 2  dedup= 3  apriori= 3  axis= 3   rho(dedup,apriori,axis)=0.966,0.986,0.965
  y_quant_uninformative_flex192_c10          causal  63/210  all= 2  dedup= 3  apriori= 3  axis= 3   rho(dedup,apriori,axis)=0.967,0.986,0.966

PER-GENE RATIO p_rule / p_all, all arms pooled: quantiles
  dedup    0.0724  0.685  0.956  1.02  1.08  1.45  5.18
  apriori  0.18  0.761  0.964  1.01  1.05  1.31  5.18
  axis     0.0724  0.681  0.954  1.02  1.08  1.45  5.18
```

**Reading of the run.** The control holds (max |dlog10 p| 4.75e-4 against the
six printed digits). Median 69 value-identical duplicate sets per region at
N = 50,000 (48 a-priori: the three rungs at or above the bound of 100, times
16 score transforms; 45 where a score transform has no admitted marker). No
rule moves the null ECDF outside the interval any more than `all` does (the
`y_null` 0.08-0.10 at alpha 0.05 is the already-open Cauchy-row inflation and
is present under every rule), and Bonferroni detections among causal genes
are identical or within one. Rank agreement with the shipped rule is
0.95-0.99. The per-gene ratio `p_rule / p_all` runs from 0.07 to 5.2 with the
quartiles at 0.96-1.08: the choice is invisible in aggregate and up to 5x on
a single gene, in either direction. So the decision is about DEFINITION, not
calibration or power: `apriori` is the paper's statistic evaluated only over
the rungs that can differ at this N, is computable at spec time, and cuts
25% of the sets; `all` is the paper's arithmetic verbatim and is what the
goldens pin. The row's candidate (a), an N-aware warning, costs nothing and
is compatible with either.
