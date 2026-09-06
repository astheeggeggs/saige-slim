# Handoff: the region path -- measured, de-churned, split, threaded (2026-09-06)

One session, four commits on `strip` after `7a9efcb8`, all gated (goldens all
PASS, standard tier FAIL 0 in both compare modes, `sim50k` outputs
`cmp`-identical, TSAN clean + plant for the threaded one). Plan:
`~/.claude/plans/can-you-summarise-what-shimmering-snail.md` (approved; steps
1-3 done, 4-5 not started). Numbers: `harness/bench/NOTES.md`, the last four
sections. Read those before this.

## What landed

1. **Step 1, `9cf9559a` -- the region C++ bucketed at N = 50,000.**
   `SAIGE_PROFILE_REGION` (the region twin of `SAIGE_PROFILE_MARKER`; serial
   arm only). The prediction was wrong in one place and it mattered: the
   VarMat product is 1.6% of the region C++; 42% was P1Mat/P2Mat churn --
   `resize` + zero-fill, strided row fills, truncation copies, and a
   write/read/delete round trip through `.bin` chunk files that EVERY region
   took because the ultra-rare collapse appends a chunk (the CLAUDE.md trap).
   Same commit: the `qfc` audit. SKAT's Davies routine is 17 file-scope
   statics and one global `jmp_buf`, read from the shipped `SKAT.so`. SKAT-O
   can never run in-process on threads; the fork pool is the only correct
   shape. `docs/state/region-perf-handoff.md` item 1 records it.
2. **Step 2a, `ee4f9721` -- chunks in memory.** `RegionChunkStore` keeps the
   `(P1 rows, P2 cols)` chunks and assembles VarMat from the SAME sub-block
   products the disk path issued (byte-identical by construction); spills to
   the old files above `SAIGE_REGION_CHUNK_MEMORY_MB` (default 8192, the
   tests force 0). The scratch is a static sized once per run and never
   zeroed; R passes 1 x 1 placeholders. Region C++ -29% / -38% (binary /
   quant SKAT-O), wall -16% / -21%, user CPU -8% / -12% -- the user-CPU
   prediction (~20%) over-counted, much of the churn was sys time.
3. **Step 2, `1852dfee` -- pass A split** into `prepareAndTest(i, k, scratch,
   reader)` and a serial `accumulate(i, k)` in marker order, in blocks of
   `markers_per_chunk`; `isSPAConvergeVec` a byte vector. Byte-neutral; ~2 s
   on a 160 s run for the block-column copies.
4. **Step 3 (the commit after `1852dfee`) -- `--nThreads` threads the region
   marker loop.** Gate: `test-step2-region-threads.R` 15/15, `sim50k` 4-thread
   outputs identical to 1-thread (24/24 + 9/9 after the race fix), goldens all
   PASS, tier FAIL 0 both modes, TSAN clean 0/0 at 4/8 threads, plant 95/109. Phase 1 on a `parallelFor`
   over the marker path's `MarkerStatePool`; ER deferred to the main thread;
   errors and short reads acted on in serial order; the count printed once
   from R; the SKAT-O pool started before the first C++ call; a `region`
   workload and plant anchor in `harness/tsan.sh`;
   `tests/testthat/test-step2-region-threads.R` (P2 escalation, P3
   conditional, P5 BURDEN at 8 threads, P36/P37 quantitative, all byte-
   identical with the thread count asserted as printed). Wall at 4 threads:
   SKAT-O 1.69x / 1.93x, BURDEN 1.57x / 1.64x.

## The one number that did not meet its target, and why

BURDEN was projected at 3.36x (`region-perf-handoff.md`, "93.7% C++") and
measured 1.57x. Located from the step-1 buckets: the C++ that threads is
prepare + per-marker test (~55% of the run); the ultra-rare pass (pass B) and
the BURDEN emitter (pass D) are 15 serial per-cell `getMarkerPval` calls
each, plus accumulate and the R side, ~45% serial. Amdahl with the 1.39x
hardware inflation floor predicts 1.55x. **Step 3b (not started): thread the
per-cell tests of passes B and D** -- the cells are independent until their
outputs are written; the residualisation and `genoSumMat` updates that
precede them are order-dependent and stay serial.

## Traps met

- **`run_arm.sh` names logs and outputs by (workload, trait, rep tag, nt)**;
  a BURDEN and a SKAT-O arm with the same tag overwrite each other. Use
  distinct tags per kind (the `region_*` drivers in the scratchpad do).
- **zsh does not word-split unquoted variables**: a golden command pasted into
  a zsh heredoc arrives as one argument. Run such snippets under `bash`.
- **The bundled SKAT-O quantitative model has only `null` variance-ratio
  rows**: P36/P37 run WITHOUT `--sparseGRMFile`. The first version of the
  threads test passed it and died in `Get_Variance_Ratio`.
- **A frame variable written from every worker is a race even when every
  write is the same value.** TSAN found exactly one: `hasVarRatio`, the
  out-flag of `selectVarianceRatio`. Per call now, as in `testOneModel`.
- **`sample` is the instrument for "where inside a bucket"**: the buckets
  said "37% ultra-rare pass, tests a fifth of it"; the tree said `resize` and
  `bzero`. Reading the code did not find it.
- VS Code's C++ indexer re-indexes `Main.cpp` after every edit at ~75% of a
  core; absolute seconds in the same-hour reps carry that.

## Not done

- **`docs/state/next.md`'s pointer line is NOT updated**: another session has
  an uncommitted section in that file ("REFUTED: port the repo to a faster
  language") and the shared-checkout rule is not to touch it. The owner
  should point `next.md` at this file and fold in the step-3b item.
- Step 4 (several null models per region) and step 5 (batch the SKAT-O cells)
  are as the plan describes them. Step 1's numbers size step 4: multi-model
  regions save the prepare bucket, 18-23% of the SKAT-O C++ (~10% of a
  SKAT-O run per extra model) but 35-59% of BURDEN's; the I/O pass is the
  rest of the case.
- The remaining serial churn for step 2's successor: the truncation copy into
  the chunk store (3 s) and the strided `P1Mat.row` fills + per-cell O(N)
  bookkeeping in the ultra-rare pass (9-12 s at N = 50k).
