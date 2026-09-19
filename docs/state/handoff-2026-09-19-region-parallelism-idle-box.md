# Region parallelism, measured on an idle box (2026-09-19)

Supersedes the numbers in `handoff-2026-09-19-region-threading-measured.md`,
which were taken on a 4-core laptop under load and were CPU-starved. **The
conclusion in that file -- and the prediction I replaced it with -- were both
wrong.** This is the measurement.

Box: 10 physical cores, 1-min load 2.39, **7.6 free**. sim50k BURDEN
(`--r.corr=1`), `--is_fastTest=FALSE`, PGEN, 3 reps.
`harness/skato_pvalue/regthread/bench.sh`.

| arm | wall | user | cores busy | wall speedup | CPU inflation | efficiency |
|---|---|---|---|---|---|---|
| serial | 11.06 s | 10.15 s | 0.92 | 1.00x | 1.00x | 1.00 |
| one process, `--nThreads=4` | 7.84 s | 12.04 s | **1.54** | 1.41x | 1.19x | 1.19 |
| four processes, quarters | 5.90 s | 18.74 s | **3.18** | **1.87x** | 1.85x | 1.02 |

All three reps reported `rows 3345 ... IDENTICAL`, so the split is a legitimate
parallelisation and the timings are valid.

## What was predicted, and why it was wrong

I predicted threads ~2.7x and processes ~2.1x -- processes LOSING -- by
assuming `wall = user / 4` on an idle box. Measured: 1.41x and 1.87x, processes
winning.

The assumption was false. **The threads arm reached only 1.54 cores of the 7.6
free.** It cannot use four. That is the ~45% serial fraction
`harness/bench/NOTES.md:1409-1419` located and modelled
(`1/(0.45 + 0.55 x 1.39/4) = 1.55x` predicted against 1.57x measured there):
threading markers WITHIN a region cannot parallelise the ultra-rare pass, the
BURDEN emitter, the accumulate or the R side. Four processes run four copies of
that serial section concurrently, and reach 3.18 cores.

So `next.md`'s original instinct -- thread across REGIONS, because that is what
attacks the serial fraction -- was RIGHT, and two arguments I made against it
were wrong: "processes already give the speed" (they give MORE, which is the
point) and "the loss is CPU inflation so more parallelism will not help" (the
loss is the serial fraction).

## What it decides

1. **For wall clock today, split the group file across jobs.** 1.87x at 4-way,
   measured, no code. The four quarters reconstruct the whole run exactly.
2. **Nothing reaches the ~2.9x that `next.md` advertised** -- that figure had no
   source and is most likely `4 / 1.39`. Best measured is 1.87x.
3. **The case for in-process region threading is efficiency, not peak speed.**
   Processes buy 1.87x at 1.85x the CPU and 4x the RSS (efficiency 1.02 -- close
   to pure resource purchase). Threads convert CPU to wall far better
   (efficiency 1.19) but are capped at 1.41x by the serial fraction. Threading
   the REGION loop in process would aim at 1.87x with threads' efficiency:
   ~1.2x CPU, 1x memory. **That is worth building if CPU-hours or memory bind
   -- a cluster that bills CPU, or one fat node running many cohorts. It is not
   worth building for wall clock alone, because job-splitting already has that.**
4. **The `qfc` thread-safety blocker does NOT apply here.** BURDEN never calls
   SKAT-O (`--r.corr=1`), so the file-scope statics that force SKAT-O onto a
   fork pool are irrelevant to this path. In-process region threading for
   BURDEN is feasible.

## Caveat

10-core Apple Silicon against the 4-core Intel box the earlier numbers came
from, so absolute times are not comparable across the two handoffs. The three
arms against EACH OTHER are what this measures, and that comparison is sound.
