# Prediction: the FlexRV RAM/CPU levers on the step-2 region path

Written 2026-09-12 04:50, before any gate ran, against HEAD `48ea91a8` (clean
tree; `SAIGE.so` `fd029535` installed from a `git archive` of it into a
session scratch library, the reference for every paired run below).

Question (owner's, 2026-09-12): apply the five levers the threading/RAM
assessment named, one gated commit each, comparing to the existing code.
Fixture for the 50k arm: `SAIGE-sim/harness/work/sim50k`, `y_null`
(quantitative, frozen model), pedigree-truth sparse GRM, ONE cell
`high:medium` at max MAF 0.001, BURDEN, `--flexRV_score=SIM` (192 sets),
the `harness/work/flexrv50k/run.sh` recipe. HEAD on this box today, load
3-5: 1 thread 53-60 s real / 46-47 s user / 875-877 MB; 4 threads 44-47 s
real / 68-71 s user / 877-888 MB; C++ profile at 1 thread prepare 4.7,
urPass 14.5 (serial, `urTest` 0), emit 21.0 of which emitTest 15.4
(threaded), urCollapse 0.14, total 44.9 s.

## Stage 1: levers 3 + 4 (records point at the columns; zero list derived per worker)

- Every committed golden byte-identical in both compare modes. Nothing
  numeric changes: the test reads the same doubles through a view, and
  `complementIndices` is `arma::find(v == 0)` by construction (increasing
  order, complement of the nonzero list that `find(v != 0)` /
  `imputeGenoAndFlip` produced).
- 50k output byte-identical to HEAD's at 1 and 4 threads.
- Peak RSS at K = 192 falls by the four N-length buffers per set the
  records held (2 x 400 KB x 2 kinds x 191 sets ~ 300 MB): 876 MB -> ~570 MB.
  K = 1 unchanged (~505 MB).
- CPU: `urPass` and `emit - emitTest` drop a little (one memcpy and one
  `find` fewer per set), well inside this box's noise; no change predicted
  in emitTest.

## Stage 2: lever 2 (per-set prepares on the pool)

- Byte-identical everywhere: each set's column has its own accumulators,
  the order-dependent state (`urTested`, `i1InChunk`, `burdenTested`) is
  assigned in a serial pre-pass in cell order.
- 4-thread wall on the 50k arm improves from ~1.2-1.35x to >= 2x over
  1 thread (the parallel fraction goes from 34% to ~80%).

## Stage 3: lever 1 (carrier-space collapse)

- Byte-identical: the per-(sample, set) accumulation order is the marker
  order, unchanged; `mean`/`accu` reproduced with armadillo's two-accumulator
  parity split over the carriers in increasing index order (to be verified
  against the installed armadillo's `arrayops::accumulate` before the gate).
- Peak RSS at K = 192 within ~30 MB of K = 1; urPass and D1 fall to the
  cost of the tests alone.

## Stage 4: lever 5 (twin memo on the burden test)

- Byte-identical: a served cell receives the bitwise result of a cell with
  the bitwise-equal column. The run prints an `asked / evaluated / served`
  line for BURDEN like SKAT-O's; served is a measurement, no test asserts it.
- emitTest on the 50k arm falls by ~the twin fraction (42% of rows).

## Outcome, stage 1 (2026-09-12 05:20)

Held: every golden byte-identical in both compare modes (standard tier
4424 / 4624 pass, 0 fail); the 50k output byte-identical to HEAD's scratch
library at 1 and 4 threads (and that library's output byte-identical to the
earlier run of the installed HEAD). `urPass` 14.0 -> 12.8 s, `emit -
emitTest` 5.2 -> 2.8 s at 1 thread (the memcpy and the second `find`).
**Missed: peak RSS fell 874 -> 736 MB at 1 thread (889 -> 751 at 4), 138 MB
against the 300 MB predicted.** The four N-length buffers per set are gone
from the code, so what the other ~160 MB was resident AS is not known --
the same gap the memory line's formula already had against the 2026-09-10
measurement (7 predicted, ~4.9 measured N-doubles per set). Stage 3 removes
the dense matrices themselves; what is left after it is the measurement
that will locate this.

## Outcome, stage 2 (2026-09-12 05:50)

Held: goldens byte-identical in both compare modes (4424 / 4624 pass); the
50k output byte-identical to HEAD's at 1 and 4 threads. RSS unchanged
(731 / 743 MB), as predicted. **Timing short of the prediction and not
clean:** at 4 threads 34.1 s real against 53.4 s for its own 1-thread run
(1.56x; HEAD's pair on the same alternation 47.9 / 60.1 s, 1.26x), with the
box at load 6-10 from other sessions' R processes throughout. The ">= 2x"
is not falsified by a loaded box and not confirmed either; re-time on a
quiet box before quoting a number.

## Outcome, stage 3a (2026-09-12 06:50): carrier-space storage

Held: goldens byte-identical in both compare modes (4424 / 4624 pass); the
50k output byte-identical to HEAD's at 1 and 4 threads. **Peak RSS at
K = 192: 874 -> 522 MB at 1 thread (891 -> 541 at 4)**, i.e. within ~20 MB
of the K = 1 run (~505 MB) -- what stage 1 could not locate is gone with the
dense matrices. CPU at 1 thread 45.5 -> 36.3 s user; at 4 threads 42.4 ->
24.1 s real (HEAD -> stage 3a, same alternation, load 3-6).
**The flipped pseudo-marker path** (unreachable on sim50k) is gated by a
drill: max MAF 0.5 with the collapse cutoff above 2N on the bundled binary
model, where `GENE1 missense;lof`'s pseudo-marker has alt frequency 0.778
and is flipped; SKAT-O and BURDEN arms, main and single-variant files, all
byte-identical to HEAD's scratch library -- and a planted build that drops
the flipped fill (every non-carrier's 2) DIFFERS from HEAD on exactly those
cells in both arms, so the drill sees the path (both ways).
