# Handoff 2026-09-17: the carrier-space `.bed` read -- MEASURED AND DROPPED

**Decision (owner, 2026-09-17): not merged.** The implementation, its tests
and its gate are archived on branch `bed-sparse` (`a476e478` the code,
`d6e71e43` the first version of this file), off `strip` at `6f36e781`;
`strip` was reset to that base and this file is the record. The reason is
in "Honest reading" below: the gain is 1.1-1.2x on quantitative scans and
nothing on binary, because the bed reader was never the bottleneck the
dense PGEN reader was -- the Item 3 lane-tally pre-screen (2026-09-05) had
already removed the decode of every marker the MAC gate drops. ~230 lines
and a second code path in the reader for that were judged not worth the
surface. **Do not re-open "give bed the same shape of gain as PGEN" without
a new measurement that contradicts the table below**; the branch is there
to re-gate if one appears. The pre-change install (`c05f23fb`) snapshot at
`harness/work/bedsparse-lib/base` can be deleted.

## What was asked

Item (b) of `next.md`'s ranked list after the PGEN carrier-space read
landed (`handoff-2026-09-16-pgen-sparse-read.md`): give `.bed` the same
shape of gain, because collaborators will not all convert to PGEN. The
owner chose it over dropping bed on 2026-09-17.

## What changed (5 files, +228 / -20, plus the test file)

`PlinkClass::getOneMarkerSparse` (`src/PLINK.cpp`), reached through a
`PlinkReader` override of `GenoReader::getOneMarkerSparse`
(`src/GenoReaderRegistry.cpp`) and feeding the SAME `prepareSparseMarker`
/ `imputeSparseAndFlip` on both the marker and the region path -- no
caller changed beyond the report line's wording. A PLINK-1 record is N0
two-bit codes, four to a byte, and a rare variant's record is mostly one
byte value: 0xFF (four HOM_REF codes; alt-first) or 0x00 (four HOM_ALT;
ref-first). Two passes over the N0/4 bytes, neither decoding a genotype to
a double:

1. the lane tally `countOneMarker` already used (factored into a shared
   `laneTable()`, its loop unchanged) tallies the FILE's samples, picks the
   common homozygous code as the more frequent of HOM_REF / HOM_ALT, and
   bounds the list: more than `2 * (N0/8)` file samples differing from it
   returns Unsupported, the bound pgenlib's difference lists carry
   (`PGEN.hpp`), so the sparse arithmetic runs only in the regime the pgen
   goldens pin. Unsupported means "nothing written, read as before": the
   pre-screen and the dense read follow, unchanged. The reader never
   answers Dense.
2. a byte compare skips every byte equal to the common byte; the others
   are decoded four codes at a time and only the samples that differ from
   the common code are listed, through a file-row -> analysis-index map
   built beside `m_posSampleInPlink` (`-1` for a fam sample outside the
   model; a model in another order is sorted, as PGEN's is). The last
   partial byte is walked per sample so the padding bits are never read.

The four tallies over the ANALYSIS samples (listed by code, every other
analysis sample holding the common code) go through `summariseCounts`, the
one body the dense read and the pre-screen use, so the QC gate sees the
same doubles, ref-first's exact integer complement included. Codes are
genotype codes through the allele order's map with missing as 3, so under
ref-first the common code is 2 and `imputeSparseAndFlip` folds every
marker. Const like the reads (S5): clones still share one `PlinkClass`.
`SAIGE_PLINK_SPARSE=0` is the kill switch, read once at construction with
the same expression as PGEN's.

**Reported:** the per-chunk line now reads `Sparse reads: K of M markers
in this chunk were read in carrier space.` for either reader
(`test-step2-pgen-sparse.R` greps only its prefix). P1 reads 85 of 100
markers sparsely (the other 15 exceed the bound), the same 85 under
ref-first, and 743-747 per 1,000-marker chunk on the 22chr fixture.

**Tests:** `tests/testthat/test-step2-plink-sparse.R`, 9 tests, 114
expectations, each asserting sparse == dense (kill switch) == golden byte
for byte AND that the run reported a non-zero sparse count: P1; P7 (the
22chr fam's 10,000 samples against a 1,000-sample model -- the subset map
and the padding, which the pgen read had no fixture for); P11 ref-first;
the missing fixture under best_guess / mean / minor (P42-P44);
quantitative; `--max_MAC_for_ER=200`; the pread fallback
(`SAIGE_PLINK_NO_MMAP=1`, the same bytes through a buffer); 4 threads == 1;
and the region path (P14 with its `.singleAssoc.txt`, P5 BURDEN,
single-in-group with `--is_output_moreDetails=TRUE`, 4 region threads ==
1, the profile line's `sparse=K`).

**Gate:** standard tier at `SAIGE_COMPARE_MODE=both`: FAIL 0 | WARN 2 |
SKIP 4 | PASS 5890 (5,776 + the 114 new). No golden moved -- predicted,
since the arithmetic is the merged sparse path's on new bytes. Warnings
ratchet 92 = baseline. Ledger lint, docdrift, optiondrift clean.

## Measured (sim50k, 1 thread, `step2_single` flags: `--minMAC=20`, Firth at 0.1, fastTest; base = the `c05f23fb` install; interleaved reps; every pairing byte-identical)

| arm | markers read sparsely | base user-s | new user-s | ratio |
|---|---|---|---|---|
| full spectrum (210,515), quantitative | 181,379 (86%) | 46.0 / 40.3 | 36.1 / 36.4 | **1.2x** |
| full spectrum, binary, Firth on (1 rep each) | 181,379 | 168.6 | 163.2 | 1.03x |
| MAF < 1% bed (155,576; `plink2 --max-maf 0.01`), quantitative | 155,576 (100%) | 9.8 / 9.2 | 8.8 / 8.7 | **1.1x** |
| MAF < 1% bed, binary, Firth on | 155,576 | 24.1 / 23.4 | 23.0 / 23.1 | 1.0x |

Profile (`SAIGE_PROFILE_MARKER=1`, quantitative full spectrum, buckets
summed over the 22 chunks): prepare 37.9 s dense -> 29.1 s sparse over
210,515 markers (180 -> 138 us per marker), test 12.7 -> 11.9 s over the
70,709 tested. **The PGEN reader on the same box and file: 35.5 user-s**,
output byte-identical to both bed arms. So the bed scan is now level with
the PGEN difference-list read, which was the deliverable.

## Honest reading: smaller than PGEN's ratio, for a reason that was already in the tree

PGEN's rare-only gain was 10.4x (quant) because the dense PGEN read
decoded EVERY record, including the 90% that `--minMAC=20` then dropped;
`PgenReader::countOneMarker` returns false. The bed reader has had the
lane-tally pre-screen since Item 3 (2026-09-05): a marker that fails the
MAC gate costs one tally (~5 us) and is never decoded. On the rare-only
bed the base arm was therefore already 9.8 s, not 30.9 s, and the byte
scan can only remove the decode of the 10% that pass (15,740 x ~195 us
= 3 s) while adding a scan pass to all 155,576. On the full spectrum, what
the sparse path costs per marker is dominated by the O(listed) impute /
scatter / clear on moderate-frequency markers (up to 12,500 listed at the
bound), which is the SAME for both readers, plus the test; the reader's
own pass over 12.5 KB per marker is a small part. The bed reader was never
the bottleneck PGEN's dense reader was; what this change buys is parity
and ~1.1-1.2x, and it does not make bed faster than PGEN on anything.

Two things NOT tried, each bounded small: a single pass with the common
byte guessed from the allele order (saves the tally, ~5 us x 210k = 1 s of
36); and skipping the second tally the pre-screen performs on an
Unsupported record (14% of markers, ~5 us each). Neither moves the
endpoint, which is set in `prepareSparseMarker`, shared with PGEN.

## Not done / next

- `next.md` item (b) is MEASURED AND DROPPED with this; the list's
  remaining items are (d) (deferred, pinned AoU measurement) and the
  region O(N) terms.
- The per-marker mmap page-touch is a floor both bed arms pay (12.5 KB per
  marker, 2.63 GB per full scan); not measured separately.
- The AoU exome runs use bed and pgen; this makes the bed arm's
  single-variant scan match the pgen arm's. Region arms: `prepareAndTest`
  takes the same sparse form, so the BURDEN 1.56x / SKAT-O 1.13x measured
  for pgen on 2026-09-16 should now hold for bed too; NOT re-measured here.
