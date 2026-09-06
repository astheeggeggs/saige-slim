# Handoff: is the serial step-2 writer the `--nThreads` ceiling? No. (2026-09-06, after the binary-profile handoff)

One session, no R/C++ change. The owner asked what was outstanding in the
multi-phenotype code, then for a plan to MEASURE the one adjacent item that the
work had made look more valuable -- formatting output rows on the workers --
before deciding whether to implement it. The measurement closed it.

## What landed (all untracked at the end of the session; the owner commits)

- `harness/bench_writer.sh` -- the driver: 1 binary, 1 quant, 3 quant and 5
  mixed models x 1/4/8 threads on `sim50k`, 2 timed reps each (arm order
  reversed on even reps) plus a `sample` profile at T = 1 and 8. Refuses a
  loaded or throttled machine; asserts the threaded and manifest arms PRINTED
  their counts (LEDGER #30 shape); `cmp`s every output against the T = 1 run
  and the P > 1 per-model files against the P = 1 runs; stamps commit, dirty
  flag and the loaded `SAIGE.so` md5 on every row.
- `harness/bench/writer_profile.py` -- reads a macOS `sample` tree, splits the
  main thread inside `mainMarkerInCPP` into writer / scan / wait / other
  (nearest categorised ancestor; scan frames under `parallelFor` are the main
  thread working as the Nth worker, not waiting) and lists the worker
  threads' top kernels. Reproduces the hand count on the 08-31 N = 1,000
  profile (682 / 3,066 = 22.2%).
- `harness/bench/writer_analyse.py` -- the table, Amdahl fits and four checks
  from `writer.tsv` + the profiles. Its C1 tolerance ignores row width, so
  the cross-trait spread (2.7 vs 2.0 us) reads FAIL; C2-C3 are the real
  instrument checks and pass. Its C4 "profile serial share" counts the main
  thread's own worker time as serial and is not meaningful; the NOTES section
  says what the correct reading is.
- `harness/bench/writer.tsv`; `harness/bench/profiles/writer_w1.*` (gitignored, local only;
  regenerate with the driver),
  `harness/work/writer/` (logs, outputs, `swap.log`, `build.log`).
- `harness/bench/NOTES.md` new section; `docs/state/next.md` item replaced
  with the resolution and the handoff pointer moved here.

## The result, one line

w = writer / main-thread `mainMarkerInCPP` time is 0.002-0.008 at 1 thread and
0.005-0.021 at 8, on every model count; ceiling 1.02x; **closed**. The
thread-scaling loss is CPU inflation (user CPU 1.3-2.7x at T = 8 for the same
work), not a serial section: at T = 8 the main thread is 91-95% busy as a
worker. T = 8 is no faster than T = 4 with a binary model on this 4-core i5.

## Traps met

- The installed `SAIGE.so` was four hours older than HEAD when the session
  began (built from the then-uncommitted binary-path work, by another of six
  open sessions). Rebuilt at HEAD before timing: same md5, so the committed
  source WAS the installed one, but that had to be shown, not assumed. Memory
  `saige-shared-checkout-timing-trap` records the check.
- `sample` attaches after pixi/R start-up and its cadence slows under 8 busy
  threads: main-thread samples were 57-90% of timed wall. Ratios within a
  profile are unaffected (writer samples invariant with T within 3%); do not
  read absolute seconds off a profile.
- Two desktop processes (VS Code's C++ indexer, a Chrome renderer, ~30% of a
  core each) ran through the binary cell and cost rep 2 up to 18% against
  rep 1. Paired arms within a rep, not absolute seconds, carry the verdict.
- Bash `run_in_background` has a 10-minute ceiling; a 50-minute bench must be
  `nohup`'d from a Monitor script that then tails its log.

## Follow-ups, same session

- `--nThreads` help (`extdata/step2_SPAtests.R`) says to use the PHYSICAL core
  count, with the measured cost of not doing so; `test-step2-cli-help.R` pins
  the sentence (red without it, green with it).
- The allocation-churn follow-up is CLOSED by a quiet-machine 4-thread run and
  profile (`writer.tsv` label `w2`; NOTES addendum): inflation 1.39x = the
  turbo/HT floor, allocator + lock 3.1% of self time, waiting 4.3%. Nothing
  worth a workspace refactor at 4 threads.

## Not done

- No commit (the owner did not ask for one). Files above are untracked.
- The multi-model manifest's advantage over serial runs shrank from 2.99x
  (Phase 3) to 1.55x because the per-model kernels got faster; that is a fact
  about the kernels, not a regression, and is recorded in NOTES.
