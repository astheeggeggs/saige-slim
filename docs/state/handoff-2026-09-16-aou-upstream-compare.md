# All of Us: slim vs upstream, speed and p-values (plan, 2026-09-16)

**Status: PLAN. Nothing here has been run.** Predictions are registered below
BEFORE any output exists; that is the point of the document. When the run
happens, results go in a new dated handoff and this one is left alone.

**Where it runs.** Inside the All of Us Researcher Workbench (TRE) only. The
data cannot leave, so every number here must be produced in the TRE and only
summary statistics come back. The repo can be cloned into the TRE from GitHub
(`astheeggeggs/saige-slim`), which is how the comparison scripts get there.
Nothing in this plan needs the pixi environment: `compare_pvalues.R` is base R.

## Why this run exists

`harness/sim/upstream_compare/RESULTS-2026-09-09.md` already compared `strip`
against upstream `main` 365deeff (**SAIGE 1.5.2**) at N = 50,000 on `sim50k`,
gene tests included:

| comparison | rows | Pearson on -log10 p | max abs dlog10 p | unpredicted |
|---|---|---|---|---|
| gene `y_binary` | 3,338 | 0.999972 | 0.2625 | 20 -> all attributed |
| gene `y_quant` | 3,338 | 0.999969 | 0.0568 | 3 -> all attributed |
| `.singleAssoc.txt` `y_binary` | 8,084 shared, 965 strip-only | 1.000000 | 0.0029 | 0 |

**The pass criterion there was "unpredicted departures = 0", not a correlation
threshold.** Keep that standard: 0.9999 with three unexplained rows is a worse
result than 0.99 with every row attributed.

Four things the AoU run adds that `sim50k` cannot:

1. **Real data.** Real relatedness, real case/control imbalance, real annotation
   files, a real gene-size distribution. `sim50k` is 210 genes on chr20.
2. **A different upstream version.** The existing study used 1.5.2. Production
   AoU is on **1.4.2**. The 1.4.2 -> 1.5.2 axis is entirely unmeasured.
3. **The LEDGER #131 revisit trigger** -- see below. It cannot be measured on
   any fixture in this repo.
4. **Whether slim's step 2 is drop-in on a 1.4.2 step-1 model.** If it is, step 2
   can be swapped without refitting, which is a migration path.

## Cohort and inputs

All of Us v8 exome, ancestry `afr`, chr20 for the main matrix and chr1 as the
extrapolation check. Three phenotypes, deliberately spanning trait type:

| phenotype | trait type |
|---|---|
| `BMI` | quantitative |
| `phe_250.2` | binary |
| `phe_153` | binary |

`[ASSUMPTION: the two phe_* are PheWAS codes for T2D and colorectal cancer.
Confirm against the phenotype file; only the binary/quantitative split matters
for the plan.]`

Step-1 models and variance ratios are the EXISTING ones, fitted by the
production 1.4.2 pipeline. **Both arms read the same model files.** That
isolates step 2 and removes step-1 differences as a confound.

Prerequisite check before anything else, cheap and blocking:
- the `.rda` contains a variable literally named `modglmm` (`ReadModel` loads it
  blind, no schema or version field);
- the `varianceRatio.txt` is the headerless 3-column `value type category` form
  that `Get_Variance_Ratio` (`R/readInGLMM.R:378`) parses. Upstream's parser
  (`main:R/readInGLMM.R:371`) takes the same shape, so a 1.4.2 file should read,
  but this has never been tried.

## Images

Mirrored into `gcr.io/weighty-elf-452116-c7` because All of Us requires
`--use-private-address`, and dsub refuses a Docker Hub image under it (the
worker VM has no external IP and can reach only Google APIs via Private Google
Access).

| arm | image | digest |
|---|---|---|
| upstream 1.4.2 | `gcr.io/weighty-elf-452116-c7/saige:1.4.2` | `sha256:678b13d6...` |
| slim 1.5.2-dev | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev` | `sha256:eb8de37e...` |
| upstream 1.5.2 | NOT YET MIRRORED -- see "third arm" below | -- |

Reference by digest, not tag, for anything whose result is quoted: a tag can be
re-pushed, and "slim vs 1.4.2" is only a meaningful claim if the two images are
identified.

## Part A -- speed

### What is already known, so the run does not re-measure it

- **`--is_fastTest=TRUE` costs +16.3% CPU on a group run and returns
  byte-identical output** (LEDGER #51, `harness/bench/NOTES.md:515`; three arms
  reconciled to a one-parameter per-pass cost model, all within the 4-6%
  replicate spread). The region path picks its variance ratio per marker from
  MAC and never reads the flag. The production wrapper passes it. Verify once on
  a paired AoU run, then drop it.
- **A SKAT-O group run is 43.9% C++ and 53.8% R**, of which `SKAT_davies`
  adaptive quadrature alone is 44.9% (`NOTES.md:575`, Rprof, two reps agreeing
  within 0.3 pp, with a `--r.corr=1` negative control correctly reporting 0.0%).
  So threading the C++ region loop alone is bounded at **1.49x at 4 threads** for
  SKAT-O and 3.36x for BURDEN. Slim also threads the SKAT-O post-processing in a
  forked pool, which attacks the other half; **end-to-end region thread scaling
  for slim at biobank N is NOT measured. That is a gap, not a known.**
- **bed vs pgen decode**: at N = 50,000, warm cache, 1 thread, upstream is
  faster on pgen by 9.5% while **strip is faster on bed by 9.8%** -- our PLINK
  path got 1.37x faster and the PGEN path only 1.13x, which flipped the ranking.
  Absolute numbers predate `b11e6b07`; the ranking stands.
  **But decode is not the cost function here.** dsub localizes the whole file
  from GCS before the script runs, and the `.pgen` is ~6.4x smaller than the
  `.bed`. Expect pgen to win on TASK WALL TIME via localization while losing
  slightly on decode, and expect that win to shrink as gene-test time grows.
  Upstream 1.4.2 refuses pgen without `--AlleleOrder=ref-first`; slim's reader
  does not need it. Do not let that flag difference become part of the
  comparison silently.

### Registered prediction, chr20

~0.77 s CPU per region pass at N = 50,000 on a 5-anno x 3-maxMAF grid
(`NOTES.md`). The AoU grid is 6 anno x 3 maxMAF, so ~0.92 s/gene; chr20 carries
roughly 540 genes; escalation adds ~10%. **Predicted gene-testing cost: 8-10 CPU
minutes**, plus startup and localization. If a chr20 task costs materially more
than ~15 CPU minutes, something other than gene testing dominates and that is
what to profile -- do not tune the score test, which is only 8.3% of a group run.

### Arms (chr20, all three phenotypes, 2 reps, order swapped between blocks)

| # | image | input | threads | answers |
|---|---|---|---|---|
| 1 | 1.4.2 | bed | 1 | baseline |
| 2 | slim | bed | 1 | version effect, like-for-like |
| 3 | slim | pgen | 1 | reader + localization effect |
| 4 | slim | pgen | 4 | thread scaling at biobank N |

Plus a **truncation series** on arm 2: the same group file cut to K genes,
K in {25, 50, 100, all}. Regress wall time on K; intercept = fixed cost (image
pull, localization, R startup, GRM and model load), slope = per-gene cost.
The region path has **no per-chunk `proc.time()` instrumentation** -- only the
marker path does (`R/SAIGE_SPATest_Marker.R:156`) -- so this regression is the
only way to split startup from work without a code change. Non-linearity is
itself the finding.

Then **one chr1 task as a positive control**: predict its wall time from the
chr20 regression by gene count, write the prediction down, run it once. Landing
within the replicate spread validates the two-parameter model and licenses
extrapolation to all 22 chromosomes without running them. Missing it most likely
means gene SIZE matters as well as gene count, which chr20's distribution may
not represent.

### Instrumentation the wrapper needs

GNU `time -v` is often absent from these images. Poll instead, around the
`"${cmd[@]}"` call:

```bash
"${cmd[@]}" &
SAIGE_PID=$!
PEAK=0
while kill -0 "${SAIGE_PID}" 2>/dev/null; do
  HWM=$(awk '/VmHWM/{print $2}' "/proc/${SAIGE_PID}/status" 2>/dev/null || echo 0)
  (( HWM > PEAK )) && PEAK=${HWM}
  sleep 5
done
wait "${SAIGE_PID}"
echo "PEAK_RSS_KB=${PEAK}"
```

Plus `date -Is` landmarks either side of the Rscript call, so the log carries its
own wall-clock split independent of Batch's accounting.

### Node sizing and threads vs instances

There are 66 naturally independent tasks (22 chr x 3 phenotypes). Independent
tasks scale linearly; threads hit Amdahl at ~1.5x on SKAT-O. **So 1 core per
task x many tasks is cost-optimal UNLESS memory forces otherwise**, and that is
the entire question:

- peak RSS per instance small (say <= 4 GB): pack N single-threaded instances on
  an N-core `n2-standard` and never thread;
- peak RSS large (16 GB+, plausible with a 50k sparse GRM): a 4-core node fits
  one instance, three cores idle, and threading is then worth its Amdahl loss
  only because the cores are already paid for.

**Peak RSS is therefore the first measurement, not the last.** With it, cost per
task = (machine $/hr / instances that fit) x task hours, compared across
`n2-standard` (4 GB/vCPU), `n2-highmem` (8 GB/vCPU) and `n2d` (AMD, typically
~10% cheaper).

### Pre-emptible

Viable, probably a large win, but **the restart machinery will not save you**.
Slim's region path has chunk-level restart (#140: the index line records both
output files' sizes and a restart truncates to them), but under dsub outputs are
delocalized only at task end, so a preempted task restarts on a fresh VM with an
empty working directory. A preemption costs the whole task.

So it is a straight expected-value question: Spot is 60-91% off, which pays
unless preemption probability over the task duration is high. At 10-30 minute
tasks it almost certainly pays. Use `--preemptible` with `--retries 2` and
**measure the actual preemption rate on the batch itself** (`dstat --status '*'
--full` shows retried attempts). Treat any chromosome whose task runs for hours
separately.

## Part B -- p-values

### Flags that must be pinned, or version is confounded with default

| flag | slim default | upstream default | action |
|---|---|---|---|
| `--minGroupMAC_in_BurdenTest` | 1 (LEDGER #126) | 5 | set explicitly, same on both |
| `--is_Firth_beta` | -- | -- | same on both |
| `--is_fastTest` | -- | -- | same on both |

### Controls, before any number

1. **Self-reproduction (the null).** Run one arm twice on identical inputs;
   outputs must be byte-identical. The sim50k study got 72/72. If AoU does not
   reproduce, STOP -- nothing downstream is interpretable.
2. **Negative control.** Run the comparison script on an arm against itself. It
   must report zero differences. This repo has hit the failure mode of a
   comparator that compares nothing.
3. **Predictions registered before looking at output.** That is this document.

### Registered predictions: rows that SHOULD move, slim vs upstream

| fix | expected effect | where visible |
|---|---|---|
| #28/#29 | p-values differ at AF > 0.5 markers; `AF_Allele2` identical; slim tests markers at MAC exactly `--minMAC` that upstream drops | `.singleAssoc.txt` |
| #103-#106 | burden variance changes -> gene p-values move | region file |
| #118 | `MAC`, `MAC_case`, `MAC_control`, `Number_rare`, `Number_ultra_rare` change; fractional MACs (3.34, 82.517) disappear; every weight set of a cell prints the SAME MAC as the plain unweighted run | region file |
| #130/#133 | pseudo-marker MAC no longer folded; `MAC_case + MAC_control = MAC` holds above alt frequency 0.5 | region file |
| #135 | BURDEN gating MAC is the count, not a folded frequency -- inert unless `minGroupMAC` is raised | region file |

**A statement that must hold, and is a bug if it does not:** on `--r.corr=0`
(SKAT-O) runs Firth cannot move any p-value. It is entered after `pval` is final
and writes only `t_Beta`/`t_seBeta` (`src/SAIGE_test.cpp:942-1017`), and both
branches pin `|BETA|/SE = |qnorm(p/2)|`. If a p-value moves between a Firth-on
and a Firth-off pair, that is a defect, not a finding.

### Tooling

`harness/sim/upstream_compare/compare_pvalues.R` -- base R, no packages, runs
anywhere in the TRE:

```
Rscript compare_pvalues.R gene   <upstream.txt> <slim.txt> <out.tsv> [label] [vr_up.txt vr_slim.txt]
Rscript compare_pvalues.R single <upstream.txt> <slim.txt> <out.tsv> [label]
```

It reports marker/gene SET differences BEFORE any statistic, then tags every
shared row with the predicted cause (`AF>0.5`, `VR`, `VR-mix`, `SKATO=1`,
`Cauchy`, `tiny`, `small`, `none`) and `UNPREDICTED` for anything else.
`compare_firth.R` reads columns as CHARACTER, so it compares printed bytes
rather than re-parsed floats -- keep that property.

Run it either inside a dsub task using either image (both carry R) or in an R
environment in the workbench. Keep `.markerList.txt` on every arm: without it an
unattributed row stays unattributed, because there is no other way to recover
which markers entered which cell.

### The third arm, recommended

Mirror `wzhou88/saige:1.5.2` as well. Then 1.4.2 -> 1.5.2 is upstream's own
changes and 1.5.2 -> slim is ours, and the latter connects directly to the
sim50k study. Without it every observed difference is the sum of two unrelated
things and attribution is guesswork.

## Part C -- the LEDGER #131 measurement

**This is the reason to keep the run's marker-level output even if everything
else comes back clean.**

#131's deferred half -- an ALT-major annotated site stays in its annotated bin,
weighted as rare on its folded MAF, with no warning -- has a revisit trigger
that is ONE number: **the fraction of annotated sites in a real callset whose
ALT allele is major.** It cannot be measured on anything in this repo. The 70
bundled hard-call markers are all below AF 0.5 (`--AlleleOrder=ref-first`
inverts all 70 rather than giving a mixture), and the one native exception,
`rs26` in the PGEN dosage goldens, is ALT-major by 3.4e-5 -- it exercises the
fold branch and separates nothing.

The AoU chr20 `afr` run produces it directly: `AF_Allele2` in
`.singleAssoc.txt` is the per-ALT frequency, and `.markerList.txt` says which
markers entered which annotated cell. The fraction with `AF_Allele2 > 0.5` among
annotated sites is the number. Extract it whether or not the rest of the
comparison is clean, and put it on the ledger row.

Worth recording alongside it: the distribution, not just the fraction. A trigger
met by sites at 0.500001 means something different from one met by sites at 0.8.

## What to write back

A new dated handoff with: the four speed arms and the truncation regression
(intercept, slope, R^2); the chr1 prediction and whether it landed; peak RSS;
the node/threads recommendation with its arithmetic; the observed Spot
preemption rate; the `compare_pvalues.R` tag tables with UNPREDICTED counts; and
the #131 fraction with its distribution. Then update `docs/state/next.md`'s
first line to name it.
