# Handoff 2026-09-12: the FlexRV RAM/CPU levers (step-2 region path)

Read `docs/state/prediction-flexrv-ram-levers-2026-09-12.md` (prediction and
outcome per stage) and `harness/bench/NOTES.md` "FlexRV RAM/CPU levers".

## What landed (one gated commit per lever)

1. `bb4003a0` levers 3+4: the persistent URCell/BurdenCell records point at
   their columns; the zero-index list is `complementIndices()` per worker.
2. `a1f44cdd` lever 2: the ultra-rare pass (B1) and the burden prepare (D1)
   run one cell per task on `runCellTests`; a serial post-pass in cell order
   assigns rowInChunk / urTested / i1InChunk / burdenTested / the refused
   count.
3. `d73a9981` lever 1 (storage): `CarrierStore` in `src/Main.cpp` -- the three
   per-(sample x cell) matrices as one q-double column per carrier; per-cell
   N-vectors scattered into three new `MarkerScratch` vectors for a call and
   cleared after; a flipped pseudo-marker carried as a per-cell fill.
4. `6ffa219a` lever 1 (CPU): the reductions over the carriers in sample
   order (parity-split sums, ranked for the case/control gathers; the
   un-flipped no-missing `imputeGenoAndFlip` reproduced statement by
   statement; a flipped cell keeps the dense routines).
5. `4535ee99` lever 5: identical BURDEN cells evaluated once (D1's post-pass keys a
   cell on its column bytes, its count bytes and its set's ER eligibility);
   `resetBurdenMemoInCPP` / `getBurdenMemoInCPP`, `R/BurdenMemo.R`, a
   `BURDEN cells: N asked for, M evaluated, K served from an identical cell`
   line once per run, and an assertion on P58's run in
   `tests/testthat/test-step2-flexrv.R` (bound, not the number).

## Gates (every stage)

Goldens byte-identical in both compare modes; sim50k FlexRV k192 BURDEN
byte-identical to a scratch install of `48ea91a8` at 1 and 4 threads; the
flip drill byte-identical (SKAT-O and BURDEN arms) with its plant caught.
Result: peak RSS 874 -> ~520 MB (K = 1 is ~505), wall 57 -> 19 s at 1 thread,
41 -> 13 s at 4, on this loaded box; 41% of the 50k run's burden cells served
from a twin.

## Not done / open

- The 4-thread scaling was never measured on a quiet box (load 3-10
  throughout). Re-time before quoting a number.
- The binary FlexRV arm on sim50k was not run (the flip drill is the only
  binary evidence at scale; the bundled goldens P58/P60 are binary at
  N = 1,000).
- SKAT-O at K = 192 keeps the dense per-model residualised copy of
  `genoSumMat`, the P1/P2 scratch and the chunk store (N x cells each); the
  recommended FlexRV run stays `--r.corr=1`.
- The marker loop's phase-2 accumulate is still serial; it is cheap now
  (`urCollapse` 0.04 s) and was not touched.
- `bash harness/prepush.sh` has not been run by this session; the standard
  tier passed in both modes after every commit.
- The drill scripts (`run_lib.sh`, `flipdrill.sh`, `parity.cpp`) live in the
  session scratchpad; the flip drill is worth a characterization test
  (max MAF 0.5, `--MACCutoff_to_CollapseUltraRare=2000` on the bundled
  binary model reaches the flip; no golden exists for it).
