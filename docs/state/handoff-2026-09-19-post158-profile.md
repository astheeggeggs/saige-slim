# Post-#158 profile, measured (2026-09-19)

The `next.md`/research-note item "re-profile post-#158 before anything else"
is DONE. Every ranking in `docs/state/research-2026-09-18-skato-speedups.md`
rested on a COMPOSED figure -- two separate pre-#158 profiles glued together.
This replaces it with one measurement. **The composed share was right; the
absolutes were low.**

Build `0b966bf8`, `harness/skato_pvalue/profile_region.sh`, sim50k binary,
1 thread, sequential arms on a quiet box (load 1.7; 111 s wall against 109 s
CPU = 98%, so the timing is not contended). Data
`harness/work/prof-20260919/prof/` (gitignored).

## The two arms

`--r.corr=1` (BURDEN) is the negative control: it never calls SKAT-O. Confirmed
in the output -- 0 SKAT-O cells, `varmat=0`, `flush=0.0095`, `rareOther=0.086`.

| | SKAT-O arm | BURDEN control |
|---|---|---|
| user CPU | 96.69 s | 38.51 s |
| C++ region buckets | 52.44 s | 34.68 s |
| peak RSS | 786 MB | 378 MB |
| cells | 2911 SKAT-O | 2973 BURDEN |

The control fixes the non-region overhead (startup, model load, IO) at
38.51 - 34.68 = **3.83 s**, which is what lets the R-side kernel be taken by
difference.

## Result

| component | seconds | share of user CPU |
|---|---|---|
| C++ region buckets (O(N)) | 52.44 | 54% |
| **R-side SKAT-O kernel** | **40.42** | **42%** |
| non-region overhead | 3.83 | 4% |

[ASSUMPTION: the kernel figure is `user - C++ - overhead`, not a direct
measurement. Basis: the BURDEN arm pins the overhead term, and it is the only
arm-independent quantity. A direct kernel timer would settle it.]

**#158's whole-run effect, measured: 207 s -> 96.69 s = 2.14x.** The composed
estimate was ~2.5x; it lands under, as predicted, because the shipped screen
threshold is 0.05 (81%/89% of cells skipped) and the composed tables used 0.01
(88%/96%).

SKAT-O over BURDEN costs **58.18 s (2.51x the run) and 408 MB** -- 40.4 s of R
kernel plus 17.8 s of extra C++ (`rareOther` 5.58, `flush` 3.97, `varmat` 3.37,
and the `urPass`/`rareArm` deltas). The 408 MB is the P1/P2 store.

## What it decides

The research note's ranking **survives**, and the cap on kernel work is now
measured rather than composed:

- Making the ENTIRE kernel free is worth **1.72x at N=50k**
  (96.69 -> 56.27 s), not more.
- At N=400k the O(N) buckets scale ~8.5x while the kernel does not scale at
  all (Davies never sees N): C++ ~446 s against a flat ~40 s kernel, so the
  kernel becomes **8% of the run** and the same lever is worth **1.09x**.

So items 4 and 5 of the research note (thread-safe `qfc`, the rank-2 eigen
assembly) are bounded by 1.72x at 50k falling to 1.09x at 400k, and
**carrier-space Phi keeps first place** -- it attacks the 54% that grows, and
removes the 408 MB.

## Two things this run does NOT settle

- **The `region_postprocess` double drive was not measured.** It only exists
  when the forked SKAT-O pool runs, and `canParallelSkato(1)` is FALSE, so at
  `--nThreads=1` there is no fork and no second drive. That needs a threaded
  arm; `profile_region.sh` hardcodes `--nThreads=1`.
- This fixture is NOT null: **82 of 2911 cells (2.8%) took the saddlepoint
  tail**, matching the 82/1645 recorded on 2026-09-17. Useful for the
  signal-rich cost question -- it is a real datapoint at 2.8% tail occupancy,
  against the ~20-30% the crossover model puts break-even at.
