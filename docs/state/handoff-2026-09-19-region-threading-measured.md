# Region-level parallelism measured; the ~2.9x ceiling is not reachable (2026-09-19)

`docs/state/next.md` has carried, since `ddedfa6d` (2026-09-07): "Still open,
unscheduled: if BURDEN at scale matters, thread across REGIONS rather than
markers, which approaches the hardware ceiling (~2.9x at 4 threads here)."

**That number has no source.** It appears once in the tree, landed with no
supporting measurement, and is most likely `4 / 1.39` -- perfect 4-way scaling
divided by the one hardware inflation term that WAS measured
(`harness/bench/NOTES.md:1240`). It is an arithmetic upper bound that assumes
the serial residue goes to zero, not an observation. This replaces it with a
measurement.

Prior art, checked first: region-level threading has **never been attempted** --
three forward-looking sentences, no code, no abandoned branch (`git log --all
-S`). The cross-region idea IS already half-implemented for SKAT-O: step 5
(`1c706bd2`) batches a whole batch of regions' SKAT-O cells into one pool call,
which took binary SKAT-O from 1.57x to 1.79x/1.82x. What is unparallelised
across regions is the C++ region body, which is the whole of BURDEN. So this is
a BURDEN-path question only.

## The measurement

`harness/skato_pvalue/regthread/bench.sh`. sim50k, N = 50,000, 210 genes,
BURDEN (`--r.corr=1`), `--is_fastTest=FALSE`, PGEN reader, 2 reps, quiet box
(load ~3). Three arms against one serial baseline; the third splits the group
file into four quarters (53/52/53/52 genes) and runs four one-thread processes
concurrently -- i.e. exactly what a cluster job-split does.

| arm | median wall | vs serial | peak RSS |
|---|---|---|---|
| serial, 1 thread | 29.1 s | 1.00x | 357 MB |
| one process, `--nThreads=4` | 19.8 s | **1.47x** | 356 MB |
| four processes, one quarter each | 15.7 s | **1.85x** | **1,429 MB (4.0x)** |

**Correctness control, run because the timings are meaningless without it:**
the four quarters produce 845 + 829 + 842 + 829 = 3,345 rows, exactly the
whole-file run's 3,345, and the sorted content is IDENTICAL. Regions are
genuinely independent and the split is a legitimate parallelisation. (The first
attempt of this benchmark failed outright -- I had dropped the genotype file
arguments -- and its merge check reported "IDENTICAL" while comparing two EMPTY
files. The row-count guard above exists because of that; an equality test
defaults to passing when its inputs are missing.)

## What it decides

1. **Processes beat threads, 1.85x against 1.47x** -- a further 1.26x. The
   job-split route is the better one, and it needs no code.
2. **Nothing here approaches 2.9x.** Four whole processes share no serial
   section at all and still reach only 1.85x of a possible 4x, so the residual
   loss is the hardware inflation floor and per-process startup, not a serial
   fraction that region threading could remove. The 2.9x is unreachable by ANY
   parallelism on this box.
3. **The case for in-process region threading is MEMORY, not speed.** Processes
   already deliver the speed; what they cost is 4x the RSS (1,429 MB against
   356 MB), because each carries its own model, genotype handles and buffers.
   At biobank N that multiple is the binding constraint, not wall clock.

**Recommendation: do not build in-process region threading for speed.** If
BURDEN throughput at scale matters, split the group file across jobs -- measured
1.85x at 4-way here, near-linear on a cluster where the processes do not share
a box. Revisit only if the 4x memory multiple becomes the blocker, and frame it
as a memory optimisation then.

Caveat on comparability: the 1.47x measured here for threads is below the
1.57-1.84x `NOTES.md` records for BURDEN, because this run uses
`--is_fastTest=FALSE`, `--is_single_in_groupTest=TRUE` and the PGEN reader on a
different day. Use the three arms against EACH OTHER, which is what they were
run for; do not compare them to NOTES' table.
