# Handoff: step 5 -- the SKAT-O cells batched across regions (2026-09-06)

One commit on `strip` after `2ed00165` (step 4b): `1c706bd2`. Plan:
`~/.claude/plans/can-you-summarise-what-shimmering-snail.md` (steps 1-5 now
done). Previous handoff: `docs/state/handoff-2026-09-06-region-multimodel.md`.
The brief this closes: `docs/state/region-perf-handoff.md`, item 3.

## What landed

`R/SkatoParallel.R` used to feed the fork pool once per region -- ~15 cells,
`clusterApplyLB`, then a barrier before the next region's C++ call -- which
held the pool at 1.30x against a 1.43x ceiling and made 8 workers no better
than 4. Now the cells of a whole BATCH of regions go in one call:

- `skatoCollect(args)` drives `region_postprocess` once to gather the cells
  it would test; `skatoEvaluate(cl, jobs)` answers a batch in one
  `clusterApplyLB` call, dispatched expensive-first by `ncol(Phi)`
  (longest-processing-time-first; cell cost spans 0.5-540 ms and the top
  10% of cells hold 48% of the work) and puts the answers back in job order;
  `skatoServe(args, answers)` hands them into the same drive with the count
  checked both ways. `region_postprocess_parallel` is the one-region
  composition, kept for the unit tests.
- `SAIGE.Region`'s per-region body is staged over a region state record `S`
  instead of frame variables: stage A per region (prepare, C++ fast pass,
  stash the fast-pass single-variant rows, collect the first-pass cells);
  stage B one pool call per batch; stage C per region in order (serve,
  decide escalation, C++ accurate pass for the escalated models, append the
  stashed rows for the others, collect the escalated cells); stage D one
  pool call, then per region in order serve and finish (marker list, result
  rows). A batch is a group-file chunk, or fewer regions when the held C++
  output lists pass `SAIGE_SKATO_BATCH_MB` (default 512; ~7-10 MB per region
  and model at N = 50,000, the N x cells `genoSumMat` dominating).
  `--nThreads=1` and BURDEN take stages A, C, D region by region with one
  `region_postprocess` drive and no pool.
- The fast pass writes each model's single-variant rows to a temp file that
  the NEXT fast pass truncates, and the copy used to happen at decision time
  from that file. Two C++ exports split `copy_singleInGroup`:
  `read_singleInGroupTemp(k)` right after each fast pass,
  `append_singleInGroup(k, lines)` at decision time (`getline` strips the
  newline, the append writes one back: the same bytes).
- Each batch prints `SKAT-O batch: N cell(s) from M region(s) to the pool`
  (and `... escalated cell(s) ...` when any model escalated), so a test can
  assert the batching happened rather than infer it from bytes.

## Gate (`harness/work/region_profile/gate_s5.log`)

`test-unit-skato-batch.R` 13 (the three pieces against the serial kernel on
a 2-worker fork pool and a mocked drive; the budget's env parsing);
`test-step2-skato-parallel.R` 24 (two regions in one batch and P6's bytes;
two batches of one at `SAIGE_SKATO_BATCH_MB=0` and the same bytes; P3; P35
with both regions escalating, its `singleAssoc` rows in region order);
`test-step2-region-threads.R` 38; `test-step2-region-multimodel.R` 93;
goldens all PASS; `prepush.sh` FAIL 0 both modes (2889 / 3089).

**A trap the gate found in its own helper.** `skato_run` listed
`--is_fastTest=FALSE` in its defaults and let callers append
`--is_fastTest=TRUE`; `optparse` keeps the FIRST value of a repeated option,
so the P3 "fastTest" run never was one -- and passed, because P2, P6 and P35
are byte-identical (LEDGER #51: the region path never reads the flag). The
helper now drops a default whose name appears in `extra`, and the new P35
test asserts the printed "Non-fast test is performed" count.

## The prediction and the measurement

Prediction, written in `harness/bench/NOTES.md` before any run: pool-only
1.30x -> 1.38-1.41x at 4 workers; on HEAD with the threaded C++ loop 9-11 s
off a ~121 s binary SKAT-O run at 4 threads (7-9%); 8 workers now beating 4;
`--nThreads=1` unchanged; everything byte-identical.

Measured (`harness/work/region_profile/ab_s5.log`, `harness/bench/NOTES.md`
step-5 section; HEAD scratch install vs the step-5 build, arms alternated,
load 3.5-5.9, every file `cmp`-identical, 30/30):

| trait | threads | HEAD | step 5 | delta |
|---|---|---|---|---|
| y_binary | 4 | 106.6, 105.5 s | **92.5, 93.4 s** | **-13.1 s (-12.4%)**, 1.57x -> 1.79x over 1 thread |
| y_binary | 8 | 104.3 s | 91.5 s | -12.8 s; 8 still no better than 4 |
| y_binary | 1 | 166.4 s | 166.2 s | 0.0%: the negative control |
| y_quant | 4 | 75.9 s | **63.2 s** | **-12.7 s (-16.7%)** |

The gain is a little above the prediction (13 s against 9-11 s), user CPU is
unchanged within 1%, the 1-thread control is exact. The one prediction that
failed: 8 workers do not beat 4, so the remaining distance to the ceiling is
not the per-region barrier (see NOTES.md for the candidates).

## chr20, the middle-tier fixture

Built this session in one clean run (`CONFIG=chr20 bash harness/sim/run_all.sh`,
75 min; `harness/work/sim`, 6.5 GB, git-ignored; facts in
`harness/bench/NOTES.md`): 10,000 samples, the whole chr20, 550 genes, five
traits, 21-47 fast-pass escalations per trait. It is the fast gate for
multi-model region work from now on -- seconds on the bundled fixture
cannot reach a 100-region chunk with escalations, an hour on sim50k can.

`region_s5_chr20.sh` at 4 threads (`harness/work/region_profile/chr20_s5_nt4.log`):
the five-trait manifest against the five single runs, SKAT-O and BURDEN,
**30 / 30 files `cmp`-identical** (region file, `singleAssoc`, marker list per
model); every SKAT-O run reports 12 batch lines (six 100-region chunks, a
first-pass and an escalated call each) and 21-47 escalations per single run,
140 in the manifest. Wall at 4 threads, load 5-8 (not a timing): single
SKAT-O runs 123-184 s, the manifest 515 s (max RSS 1.5 GB); BURDEN 17-20 s
single, 30 s the manifest.

## State for the next session (written for a context clear)

- Commits on `strip`: `1c706bd2` (step 5), then the handoff commit.
  `docs/state/next.md`'s pointer line names this file; the rest of that file
  has another session's uncommitted section -- do not commit it wholesale
  (`harness/work/region_profile/drivers/stage_next_pointer.sh` stages the
  pointer line alone).
- The installed library is the step-5 build. A scratch install of `2ed00165`
  (unpatched HEAD before step 5) is at
  `/private/tmp/claude-501/-Users-dpalmer-Repositories-SAIGE/b16c097a-c66f-42e3-96d9-af2b303d045a/scratchpad/s5head/lib`
  (session-local; rebuild per `~/.claude` memory `saige-drill-library-setup`).
- Drivers of the day, all under `harness/work/region_profile/drivers/`:
  `region_s5_gate.sh` (the gate above), `region_s5_chr20.sh` (five-model
  manifest vs five single runs on the chr20 fixture, SKAT-O and BURDEN, cmp
  per file), `region_s5_ab.sh` (paired A/B on sim50k).

## Not done / caveats

- Jobs are still one cell each; the per-task master overhead (8-20 ms per
  region measured for the pool) is what remains between the measured figure
  and the ceiling. Jobs as LISTS of cells is the next lever if it is wanted.
- The batch holds every region's C++ output list until its answers are back:
  memory rises by up to `SAIGE_SKATO_BATCH_MB` per run (five-model manifests
  hit the budget sooner, so their batches are shorter).
- The SKAT-failure sidecar (LEDGER #57) records in serve order, which is
  region order within a batch: unchanged.
