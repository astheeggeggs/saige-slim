# Handoff 2026-10-03: LEDGER #176 FIXED -- the stretch gate; what the AoU sample measured; what is next

Read `python3 harness/tools/ledger.py 176` (its STATUS is the authoritative history, oldest first),
then `docs/state/handoff-2026-10-02-176-review.md` (the review that replaced the first handoff's
repair). This document is the record of the two builds and the measurement between them.

Counts here are counts of CELLS, JOBS or GENES. No carrier count, no per-carrier quantity and
no participant-level value is written, here or in any sidecar the builds produce.

## The two commits, both on `slim/strip`

| commit | what | image |
|---|---|---|
| `cbb9d6d7` | the stretch gate as an INSTRUMENT: evaluated, counted, named; no p-value changed | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-cbb9d6d7` (IMG11 in the workbench) |
| `152ffd8c` | the gate DECIDES: a cell at or above 0.5 reports the exact convolution at the unstretched score | `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-152ffd8c` (pushed 2026-10-03, gcr only; IMG12) |

Both went through the pre-push gate in both compare modes (7,058 and 7,073 passes). The first
moved no result golden over two full regenerations; the second moved exactly P69 (six cells),
P78 (two) and P77 (three pseudo-markers and the SKAT-O cells they feed), each value predicted
from the unpatched build's `SAIGE_SPA_DUMP_BELOW=1` dumps through `harness/tools/spa_exact_ref.py`
before the C++ existed, and reproduced to the printed digits.

## What the gate is

`q = Tstat/sqrt(varRatio) + m1` is compared against a CGF built on the unscaled vector, so the
variance ratio stretches the observed score by `|Tstat| |1/sqrt(VR) - 1|`. The gate is that
stretch over the sd of everything but the heaviest carrier's term, evaluated for every weighted
set's binary cell whose test reaches the saddlepoint (`TestContext::stretchGate`, both region
paths). Below 0.5 nothing changes. At or above it the cell reports the exact two-sided tail at
the UNSTRETCHED score (`src/ExactConv.hpp`: a convolution over the carriers on a lattice with a
Gaussian part for everyone else -- the fast kernel's own model with the tail computed instead
of approximated; no relatedness correction, as for a routed cell).

Switches, each read once per process: `SAIGE_STRETCH_GATE=0` (kill switch, pre-instrument
output), `SAIGE_STRETCH_GATE_APPLY=0` (instrument only), `SAIGE_STRETCH_GATE_THRESHOLD`,
`SAIGE_STRETCH_GATE_FLOOR` (listing floor, 0.1), `SAIGE_EXACT_CONV_SHADOW=1` (the convolution
with its bracket for EVERY evaluated cell; slow). Every region run prints
`Stretch gate: N weighted binary cell test(s) reached the saddlepoint; K at or above 0.5 (...)
K of them report the exact convolution ...` and names cells at or above 0.1 in
`<out>.stretchGate.txt` with `gate`, `z`, `p_spa`, `p_exact` and `reported`.

## The measurement (instrument image, seven unfiltered chromosome runs)

`fxrate2-<anc>-<trait>-c<chr>` jobs: amr EFRMB chr19 (shadow mode) and amr/eur x EFRMB, Asth,
CAD on chr20 (no shadow). Outputs in `<anc>/flexrv_caprate2/`; the counts come from the logs.

| | cells | share |
|---|---|---|
| burden cells | 732,374 | |
| reached the saddlepoint | 58,713 | 8.0% |
| gate at or above 0.1 | 162 | 0.022% |
| gate at or above 0.5 | 67 | 0.009% |
| gate at or above 1 | 4 | 0.0005% |
| gate at or above 3 | 4 | the known gene, `ENSG00000213030` |

- Nothing between 1 and 3 on any run: the damage is absorbed or catastrophic.
- The 63 cells at 0.5-1 are 28 in amr Asth chr20 (one gene's weight sets, by the look of it),
  2 in amr EFRMB chr20, 2 in eur EFRMB chr20 and 31 on chr19. Their exact p is not known: the
  shadow run's sidecar was not delivered (see the wrapper note below).
- The shadow kernel cost 2,290 s on chr19 (~140 ms per evaluated cell): on every evaluated cell
  it would double or triple a job; on the gated cells it is seconds a chromosome.
- Base durations on the previous image, for scale: amr chr20 5-8 min, eur chr20 15-23 min,
  amr chr19 13 min.

## The fixtures (variance ratio 0.942: a 3 % stretch, three times AoU's)

- P78: 4 cells reach the saddlepoint, 2 are gated (0.840, 0.765): 0.0311 -> 0.0540, 0.0299 -> 0.0527.
- P77: the ultra-rare pass, 3 gated pseudo-markers (1.08, 1.02, 1.02): 0.0438 -> 0.0465 etc.
- P69 (the FlexRV atom fixture): 109 cells reach the saddlepoint, 6 gated -- the five `sbeta5`
  sets at 1.6-1.7 had printed 0.0019-0.0050 in a COMMITTED golden where the exact tail is
  0.0255-0.0343 (7 to 14x too small), and `sbeta3_mafbeta5` at 0.57 moves 2 %. Below the gate the
  saddlepoint is ~0.6x of the exact tail on cells whose heaviest term fired: the method's own
  error, not the stretch's, left as is.
- The artefact in miniature on the atom fixture now needs BOTH nets off
  (`SAIGE_EXACT_BY_WEIGHT=0 SAIGE_STRETCH_GATE_APPLY=0`); `test-step2-exact-by-weight.R` says so.

## Decided, and not done

1. **No relatedness correction on gated cells**, like routed cells. The p-scale form
   (`p' = 2 Phibar(Phibar^-1(p/2)/sqrt(r))`) gives 0.93-0.94 of the unstretched answer on every
   synthetic shape and would move about a third of sas burden cells; it is its own decision and
   should get its own ledger row.
2. **The saddlepoint's own ~0.6x error** on cells whose heaviest term fired (seen in the real
   dumpcheck table and on P69) is the method's, not #176's. Its own row too, if anyone wants it.
3. **The R-side rescale** (`get_newPhi_scaleFactor_traitType`) keeps the same stretch on the raw
   burden vector; its effect is bounded by the `r = min(VarQ/VarQ_2, 1)` cap
   (`harness/tools/spa_176_rside_rescale.R`). SKAT-O arm and conditioning markers only.

## Next

1. ~~Push the fix image~~ -- done 2026-10-03 (`1.5.2-dev-152ffd8c` on gcr; Docker Hub not yet, it is
   not a production candidate until the binary-arm pass is read).
2. **One pass of the binary FlexRV arm on the fixed image** -- STARTED 2026-10-03 and restarted:
   `submit_v8` is `submit_v7` with IMG12, `fxv8-` names and `flexrv_v8/`; the first chunk
   (`fxv8-amr-efrmb-b0`, 5 chromosomes, on-demand) landed and passed its checks (51,530 cells
   reached the saddlepoint, 20 gated and all 20 reported exact, sidecar delivered); the rest of
   the first launch was cancelled and is to be relaunched on SPOT (`spot=True`), all five
   ancestries, ~590 jobs as v7. The splice control accepts gated cells since `8bd0d358`
   (`.stretchGate.txt.gz`, `reported == exact`); the workbench copy of `aou_spa_splice.sh` must
   be that version (its self-test says "9 positive controls ... one via the stretchGate
   sidecar"). `submit_check_v8`: old side `flexrv_v7_checked` (else `flexrv_v4_fixed`, else v4,
   by looking), new `flexrv_v8`, out `flexrv_v8_checked`; read `cells reported by the exact
   convolution at the stretch gate (LEDGER #176): N` and the absence of `ERROR`.
   **The pass ran 2026-10-03/04 (Spot, 589 jobs; three eas jobs that the first cancel had left in
   the registry were resubmitted on demand).** Gate totals from the logs, four ancestries in:
   sas 7 gated / 0 at or above 1 / 0 at or above 3; afr 974 / 192 / 5; amr 2,162 / 752 / 73;
   eur 789 / 152 / 2 -- about 3,900 cells reported by the exact convolution, 80 in the class
   where the old p was wrong by orders, most of them amr; 'reported exact' equalled the count
   at or above 0.5 everywhere and every trait ran 22 chromosomes. So the sample's 'nothing
   between 1 and 3' was a property of seven chromosomes, not of the arm. `check8-` jobs
   submitted for the four; eas to follow.
   eas: 36 gated / 0 / 0, so the arm has 3,968 gated tests, 1,096 at or above 1, 80 at or above 3.
   **The `check8-` control (2026-10-04)**: 589 checks, every one OK after two rounds -- the first
   round refused 55 (25 afr, 26 amr, 4 eur) because the wrapper's concatenated sidecar can carry
   a 9-column header over 11-field rows and the control looked `reported` up by name only
   (fixed in `aou_spa_splice.sh`, by position when the header lacks it; planted as `sghdr`).
   Accepted gated keys (Region|Weight, so fewer than tests): sas 4, eas 19, afr 450, amr 1,025,
   eur 340. `flexrv_v8_checked/` is the binary arm's final.
   **What the sidecars say (`gated_genes`, 2026-10-04; gene identifiers with trait, cell counts
   and p-value ratios -- summary quantities, no counts of people):**
   - amr: 2,162 gated cells in 174 genes, 23 traits; median correction 0.13 orders. At or above a
     gate of 3: 73 cells in 7 genes. Orders the old p sat below the exact one: EFRMB
     ENSG00000213030 14.9 (the known gene); IFHern ENSG00000182816 7.5 (invisible to the
     discordance screen); BenIntNeo ENSG00000196109 4.2; MatHem ENSG00000186971 4.0; PAD
     ENSG00000181698 3.9; MatHem ENSG00000205423 3.6; BenCervUterNeo ENSG00000196109 2.9; MatHem
     ENSG00000183161 1.8. The traits are the ones amr had prior form with (the three #172-type
     spurious hits, the flagged MatHem model). Why amr's models concentrate steep sets more than
     eur's is NOT answered here.
   - afr: 974 gated in 91 genes, median 0.19 orders; ILDSarc ENSG00000167523 3.5 orders (one
     cell); BenIntNeo ENSG00000215454 at a gate of 3.6 but only 0.32 orders -- a high gate whose
     damage did not materialise (the heavy carrier's outcome did not put the score at the atom).
   - eur: 789 gated in 60 genes, median 0.09 orders; HF ENSG00000237289 2.9 orders.
   So about ten gene-trait results in the binary arm were wrong by three orders or more, nine of
   them amr; everything else gated moved by a factor of a few at most, in both directions.
   **Summaries and the discordance screen on v8 (2026-10-04)**: `summary_v8/` built from
   `flexrv_v8_checked` for all five (provenance counts 6/15/124/124/320 exact; the summary
   script's sidecar list was extended first, `9cf28a54`). Screen (section 7b, reading
   `summary_v8`): binary gene-tests above 10 orders -- 0 in every ancestry; above 5 -- afr 1
   (ABCA4 for AMD), eur 5 (PKD1, BRCA1, MLH1, TERT -- the same real signals as before the fix),
   amr 0 (ENSG00000213030 gone from 13.9; IFHern's gene not in the top five). So the one
   residual outlier is resolved, nothing new appeared, and real signals did not move. The
   screen now has no positive control on v8 data; the v7 run is the record that it detects.
   Original item: (the gated cells cannot be found
   from existing output), as the v7 rerun was done; the quantitative arm is untouched. The
   wrapper now delivers `.stretchGate.txt` (the `OUT_FILE_STRETCHGATE` output in `submit()` and
   the `stretchGate:` pair in the delivery loop, added 2026-10-03 in the workbench), so every
   gated cell is named with both p-values.
3. ~~**LEDGER #171**~~ -- FIXED 2026-10-04 (`on.exit` removes any created file still empty when
   the fit exits without completing; step 2's load-time refusal measured to leave no file). The
   sidecar writer now emits `p_exact` and `reported` always. Rows #177 (the p-scale variance ratio
   for exact-test cells, OPEN, the owner's) and #178 (the saddlepoint's own ~0.6x error, DECIDED,
   leave as is) record the two decisions. Original item: step 1 creates `<prefix>.rda` and `.varianceRatio.txt`
   at startup and fills them after the fit, so every refusal leaves 0-byte files. Fix: keep the
   early create (it is the writability check), register an `on.exit` that deletes any created
   file still empty when the function exits by error, set a success flag before the return.
   Test: the refusal runs in `test-step1-fit-gates.R` (`fg_cli`) assert no 0-byte file remains.
   No golden moves. Also check whether step 2's own load-time refusal leaves an empty results file.

## Notebook notes (the workbench, not this repo)

- `submit_caprate2(anc, t, chrom, image, shadow)` and `gate_rate(names)` are session cells, not
  in the generator; `_v4(anc, t)` resolves a production record (eur keys are chunked, `-b0`).
- The wrapper captures SAIGE's output per chromosome and appends it at the end, so a running
  job shows nothing but the flush heartbeat; `tee` into the dsub stdout log would give live
  progress at the cost of interleaving under `PARALLEL_CHRS > 1`.
- Fetch logs sequentially; a thread pool is the shape of the earlier egress alerts.
