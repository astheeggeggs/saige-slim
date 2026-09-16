# Handoff 2026-09-16: the carrier-space PGEN read on the step-2 marker path

Branch `sparse-read` in the worktree `/Users/dpalmer/Repositories/SAIGE-sparse`
(off `strip` at `51d1a4b6`), scratch libraries under
`/Users/dpalmer/Repositories/SAIGE-sparse-lib/{base,new}`. Nothing here has
been committed or merged; the shared checkout and its installed package were
never touched.

## What was asked

"Someone said PGEN is much smaller AND much faster than .bed; we see no speed
difference." The research half is in the session's scratchpad report and the
memory note `saige-pgen-vs-bed-speed`; the short form: the size claim is a
rare-variant effect (difference lists), the speed claim is false for this
code because both readers expanded every variant to a dense N-vector, and no
primary source documents a format-driven compute gain anywhere. This handoff
is the implementation of the one lever that was on the table.

## What was measured before touching anything

`SAIGE_PROFILE_MARKER=1`, sim50k (N = 50,000, 210,515 markers), the
`step2_single` flags (`--minMAC=20`, Firth at 0.1, fastTest), 1 thread, on a
box carrying load 3-4 from other sessions. Whole-run bucket sums over the 22
chunks:

| arm | prepare | per marker | test | per tested | prepare share |
|---|---|---|---|---|---|
| bed, binary | 53.9 s | 256 us | 125.1 s | 1770 us | 30% |
| bed, quant | 41.0 s | 195 us | 14.1 s | 199 us | 74% |
| pgen, binary | 62.6 s | 297 us | 123.0 s | 1739 us | 34% |
| pgen, quant | 62.8 s | 298 us | 13.0 s | 184 us | 83% |

All four outputs byte-identical across readers. `prepare` is read + decode +
QC + `imputeGenoAndFlip` + the two index copies: every one of them O(N) per
marker, paid by the 140k markers `--minMAC=20` then drops as much as by the
70,709 it keeps. The test buckets were already O(carriers) on the fast path
(`scoreTestFast` over `indexForNonZero`, `SPA_fast` over the non-zero
genotypes), so the prepare bucket was the whole prize.

Region path, same fixture, pgen reader, y_binary, `SAIGE_PROFILE_REGION=1`:
prepare 17.6 s of a 53.8 s BURDEN run (33%), 16.7 s of a 211.8 s SKAT-O run
(8%; the rest is R-side Davies). NOT changed here -- see "not done".

## What changed (9 files, +630 / -41)

1. **`Geno::GenoReader::getOneMarkerSparse`** (`src/GenoReader.hpp`): an
   optional read returning `Unsupported` (PLINK; the caller reads densely as
   before), `Dense` (the reader expanded the record itself and the
   `getOneMarker` outputs are filled), or `Sparse` (a `SparseMarker`: common
   genotype CODE plus the listed analysis indices, strictly increasing, and
   their codes; the counts in `MarkerMeta`).
2. **`PGEN::PgenClass::getOneMarkerSparse`** (`src/PGEN.cpp`):
   `PgrGetDifflistOrGenovec` with `max_difflist_len = 2 * (N/8)` -- the
   raregeno buffer's capacity as pgenlibr sizes it, and plink2 `--glm`'s own
   cap. A record the writer stored as a difference list no longer than that
   comes back as (common, rows, codes); rows map to analysis indices through
   the new `PgenShared::analysisIndexOfFileSample` (the inverse of
   `posSampleInPgen`, built once), rows outside the model are dropped, and
   the list is sorted into analysis order (a no-op when the model's samples
   are in file order). Everything else -- a dense-stored record, a common
   code other than 0/2, a file with a dosage track, a multiallelic file, a
   plink1-compatible-mode file, `SAIGE_PGEN_SPARSE=0` -- takes the dense path,
   and the dense expansion is `getOneMarker`'s own loop moved into
   `expandFromReadBuffer` unchanged, so it cannot differ from it.
3. **`prepareSparseMarker`** (`src/Main.cpp`) and **`imputeSparseAndFlip`**
   (`src/UTIL.cpp`): the dense half's decisions and arithmetic over the k
   listed samples. The identities each rests on are written above the
   function: a common code of 2 always flips and 0 never does (so the
   unlisted samples are exact zeros after the fold); the imputed value is
   computed from the flipped frequency as before; `.clean()` leaves a zero a
   zero; `arma::sum` is `arrayops::accumulate`'s parity-split two-accumulator
   sum, reproduced by adding each listed value into the accumulator its
   sample index selects (the CarrierStore drill of 2026-09-12); the non-zero
   index is built in increasing sample order. The carriers are scattered
   into the worker's `GVec`, which is all-zero between markers on this path
   (`MarkerScratch::gvecDense` / `sparseScattered`), and cleared after the
   models by `testOneMarker` and, idempotently, at the top of the next
   prepare -- so every dense consumer (Firth, the non-fast SPA, the PCG)
   still reads a complete N-vector.
4. **The zero-genotype index is no longer materialised on the marker path,
   for BOTH readers.** Its only consumers were its SIZE (now the identity
   `m_n - nnz` in `getMarkerPval`, exact because the two index vectors
   partition the samples on every caller) and the ER exact test, which now
   builds the complement itself when handed an empty vector (main-thread
   only, MAC <= `max_MAC_for_ER` markers). `imputeGenoAndFlip` grew a
   defaulted `t_buildIndexZero`; the region and conditional paths still
   build it.
5. **The run reports the path**: `Sparse reads: K of M markers in this chunk
   were read as pgen difference lists.` per chunk, tallied from a
   `MarkerResult::readSparse` flag by the main thread like the Firth
   counters. `tests/testthat/test-step2-pgen-sparse.R` asserts the line
   (present with K > 0 on a pgen run, absent under the kill switch).

## Gates

- `bash harness/run_tests.sh standard` at `SAIGE_COMPARE_MODE=both` against
  the patched library: **FAIL 0, PASS 5737, 521 s**. Every PGEN golden
  byte-identical with the sparse read on.
- `test-step2-pgen-sparse.R` (new, 5 tests, 50 expectations): goldens P13 and
  P20 reproduced; sparse == dense (kill switch) byte for byte under
  `best_guess`, `mean` (fractional imputed values -- the parity-split sum)
  and `minor` on the missing-genotype fixture, on the quantitative model, on
  the ER branch (`--max_MAC_for_ER=200`), and at 4 threads == 1 thread.
- Warnings ratchet: 92 src-anchored, at the baseline.
- The tier was run on the build BEFORE the idempotent top-of-prepare clear
  (a bool test per marker, unreachable unless a test throws); the rebuild
  and the tier rerun on the final source are the last step of the session
  and are recorded at the bottom of this file.

## Measured (sim50k, pgen, 1 thread, interleaved, 2 reps each, same box state)

| trait | unpatched user-s | sparse user-s | ratio | elapsed |
|---|---|---|---|---|
| quantitative | 71.4 / 70.8 | 37.3 / 38.8 | **1.87x** | 79.3 -> 43.6 |
| binary | 179.0 / 179.4 | 147.8 / 147.2 | **1.21x** | 195.8 -> 160.4 |

Outputs byte-identical base vs new and rep vs rep. 170,520 of 210,515 markers
(81%) came back as difference lists. Against the Amdahl bound from the
profile (prepare 83% / 34%) the quantitative scan realised about half the
ceiling -- the remaining prepare cost is the per-marker fixed work (the
metadata strings, the row fields, the `conv_to`) and the 19% of markers
above the cap that still expand -- and the binary scan is where the profile
said it would be, test-bound.

## Not done, in the order the numbers rank them

1. **The region path.** `mainRegionInCPP` still calls `getOneMarkerThrough`
   (`src/Main.cpp` ~2731) and `imputeGenoAndFlip` (~2747) per marker and then
   scatters into `CarrierStore`; prepare is 33% of a BURDEN run and 8% of a
   SKAT-O run on the pgen reader. The sparse read would feed CarrierStore's
   carrier list directly. Same invariants, second call site.
2. **The bed reader.** `PlinkClass::getOneMarker` walks N samples; a scan
   that skips bytes equal to the common 4-genotype pattern and decodes only
   the rest would give `.bed` the same shape of gain for rare variants
   (prepare was 256 / 195 us per marker there). `countOneMarker`'s
   byte-lookup pre-screen already exists.
3. **Dosage-carrying PGEN**: unsupported by design (the sparse form carries
   no dosage track); `PgrGetDMaybeSparse` is the API if it is ever wanted.
4. **Step 1** re-encodes each PGEN variant into the packed bed store and its
   kernels are O(N) by nature; not a lever.

## Traps for whoever picks this up

- The sparse path fires only on records the WRITER stored sparsely. A
  plink1-compatible-mode `.pgen`, or one written with a dosage track, gets
  nothing, silently except for the missing `Sparse reads:` line. The file
  size says which kind you have: sim50k's is 408 MB against a 2.63 GB bed.
- `max_difflist_len` must be below `raw_sample_ct` and within the buffers
  `initReader` carved out; `2 * (N/8)` is both. Raising it to `N/4` is NOT
  the same number when N is not a multiple of 8.
- The reader decodes ALL file samples (`_subset_size == file_sample_ct`) and
  maps to analysis order itself, so the list's ids are pgen ROWS. If a
  sample subset is ever pushed into pgenlib, the ids become subset-relative
  and `analysisIndexOfFileSample` is the wrong map.
- `MarkerScratch::GVec` is uninitialised at construction; `gvecDense` starts
  true so the first sparse read zeroes it. Do not "optimise" that away.

## Final gate (on the rebuilt library, with the top-of-prepare clear)

`bash harness/run_tests.sh standard` at `SAIGE_COMPARE_MODE=both`,
`R_LIBS=SAIGE-sparse-lib/new`, from the worktree: **FAIL 0 | WARN 2 | SKIP 4 |
PASS 5737, 507 s**; `ledger lint OK: 157 rows`, `docdrift OK: 21 rules over 39
files, 0 hits` (this handoff, the CLAUDE.md paragraph and the next.md
paragraph included), optiondrift selftest OK. The two WARNs and four SKIPs are
the tier's standing ones (Guard Malloc and the platform skips), unchanged.

## What is left in the scan, by symbol (after the change)

`sample` at 1 ms over 20 s of the patched quantitative sim50k scan (14,760
samples inside `mainMarkerInCPP`; the profiled prepare bucket on the same
run was 30.5 s = 145 us per marker, against 67.7 s = 321 us unpatched):

| symbol | samples | share of the loop |
|---|---|---|
| `PgenClass::getOneMarkerSparse` (all) | 5,221 | 35% |
| of which `expandFromReadBuffer` -- the DENSE loop for the 39,995 dense-stored records | 4,313 | 29% |
| of which `PgrGetDifflistOrGenovec` (pgenlib's read, both kinds) | 431 | 3% |
| of which `Dosage16ToDoubles` | 283 | 2% |
| `getMarkerPval` (the test) | 4,683 | 32% |
| `conv_to` (the non-zero index copy) | 178 | 1.2% |
| `imputeSparseAndFlip` | 131 | 0.9% |
| `writeOutfile_single` | 72 | 0.5% |

So the carrier-space arithmetic costs nothing measurable and pgenlib's read
is 3%; what remains of prepare is the 19% of markers whose records are dense
(non-modal count above N/8, i.e. MAF above ~6%), each walked sample by
sample through `posSampleInPgen` -- a permutation here, since the model's
sample order is the phenotype file's, not the psam's (50,000 of 50,000 rows
differ), so every dense record is a random gather over a 400 KB buffer --
with an UNCONDITIONAL `t_indexForNonZero.push_back(i)` per non-zero sample
that the caller then discards (`imputeGenoAndFlip` rebuilds the index;
`PlinkClass::getOneMarker` already dropped its copy of this push in Phase
4). Dropping that push is byte-neutral and is the next marker-path item on
the PGEN side; it needs its own gate. The permutation gather is inherent to
the model owning the sample order. The quantitative test bucket (32%) is
`scoreTestFast`'s O(k p) over those same common markers.

The sort in `getOneMarkerSparse` DOES run on this fixture (the model order
is a permutation of the file order) and is invisible in the sample, as the
k log k over tens of carriers should be.

## The question that followed: a pgen filtered to rare variants

`plink2 --pfile chr20_sim_pgen --max-maf 0.01 --make-pgen`: 155,546 of the
210,515 variants, **4.9 MB** against the full file's 408 MB (and the bed's
2.63 GB). Same `step2_single` flags (`--minMAC=20`, so 15,740 markers
tested), 1 thread, base vs patched, one run each, outputs byte-identical in
every pairing, every marker read as a difference list:

| arm | unpatched user-s | sparse user-s | ratio |
|---|---|---|---|
| quantitative | 30.9 | **3.0** | 10.4x |
| binary, Firth on (`--pCutoffforFirth=0.1`) | 43.5 | **15.7** | 2.8x |
| binary, Firth off | 33.6 | **5.5** | 6.1x |

The patched quantitative run is process startup plus ~0.5 s of scan: about
3 us per marker for read, QC, impute and (for the 10% that pass
`--minMAC=20`) the score test and the output row. With Firth off -- the
owner's plan under the Stouffer scheme (`handoff-2026-09-16-aou-first-results.md`)
-- the binary scan is within 2x of that; what separates them is the
saddlepoint on the ~5% of markers with |Z| > 2, which is O(N) per marker
(`getadjGFast`, the dot with mu). Firth on is the O(N x iterations) tail
and caps the ratio at 2.8x on this fixture.

So: **yes for single-variant scans on a rare-only pgen, an order of magnitude
for a quantitative trait and ~6x for a binary trait without Firth**, with
the caveat that these are 1-thread wall-clock-scale numbers on a fixture
where the unpatched scan was 31-44 s; at biobank N the O(N) terms grow 8x
while the O(carriers) terms do not, so the ratios should be larger there,
and NOT measured. Region tests are a different path with a different bound
(above).
