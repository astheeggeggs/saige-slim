# PGEN input for step 0, step 1 and extractNglmm (2026-09-03)

Before this, PGEN reached **step 2 only** (Phase 7). Steps 0/1 and the
`extractNglmm` diagnostic required PLINK bed/bim/fam, so plink2-native data had
to be converted back to bed just to build a GRM and fit a null model.

They now take `--plink2File` (a prefix) or explicit
`--pgenFile`/`--pvarFile`/`--psamFile`, mutually exclusive with the PLINK
options. **Every output is byte-identical to the PLINK run on the mirrored
fixture** -- step-0 `.mtx` and `.sampleIDs.txt`, step-1 `.rda` and
`.varianceRatio.txt`, across all five step-1 golden configurations.

## The design: re-encode, don't re-implement

`genoClass::setGenoObj` has always streamed the `.bed` one marker at a time into
a `ceil(N/4)`-byte PLINK-1 record and handed that to
`Get_OneSNP_Geno_atBeginning`, which does the unpack, the MAF/missingness QC,
the `round(2*altFreq)` mean imputation and the 2-bit repack.

**The pgen path produces that same record and calls that same function.**
`src/SAIGE_fitGLMM_fast.cpp` now has a `MarkerByteSource` with two
implementations -- `BedByteSource` (the original seek+read) and
`PgenByteSource` (pgenlib decode, then re-encode as a .bed record) -- and
`setGenoObjFromSource` is the single body both front-ends drive.

Byte-identity therefore holds **structurally**: there is one decode/QC/impute/
repack pipeline with two front-ends, not two pipelines that happen to agree.
The pack/unpack round trip is real overhead, knowingly traded for that
guarantee; it can be short-circuited later now that goldens pin the result.

Encoding (`Get_OneSNP_Geno_atBeginning` reads `b` = bit 2j, `a` = bit 2j+1 and
maps (b,a) -> (1,0) missing, (0,0) 2, (0,1) 1, (1,1) 0, i.e. it counts copies of
the bim column-5 allele; field value is `a*2 + b`):

| ALT count | 2-bit field |
|---|---|
| 2 | `0b00` |
| 1 | `0b10` |
| 0 | `0b11` |
| missing (NaN) | `0b01` |

## Two facts that had to be true, and were measured

1. **`modglmm` carries nothing path-dependent.** A recursive scan of
   `harness/goldens/step1/binary_sparseGRM_a.rda` finds no calls, formulas,
   closures or file paths anywhere in the object (`obj.glm.null` is a bare
   `list(model = )`). Without this, byte-identity would have been impossible on
   principle rather than on arithmetic.
2. **pvar ALT == bim column 5 for every marker.** Checked over all 128,868
   markers of `nfam_100_..._22chr` and all 10,000 of
   `plinkforGRM_1000samples_10kMarkers`, plus psam IID order == fam IID order.
   Step 1 counts copies of bim column 5; pgenlib's `alleleIdx = 1` reads ALT.
   Same quantity -- which is what makes the two readers comparable at all.
   `plink2 --bfile X --make-pgen` puts bim column 5 in ALT, so this holds for
   any fixture built that way; a pgen built with REF/ALT swapped relative to a
   bim would flip every marker. The GRM is invariant to that flip (it is a sum
   of products of standardised genotypes), so it is not a correctness hazard,
   but it WOULD break byte-identity -- which is why the check is worth
   repeating for any new fixture.

## Sample subsetting deliberately bypasses `posSampleInPgen`

`genoClass` does its own subsetting -- `ptrsubSampleInGeno` (1-based into
genotype-file row order) and `indicatorGenoSamplesWithPheno_in` (a mask over all
file rows), both built in R -- so steps 0/1 want ALL psam rows in file order.
That is exactly what `PgenClass::Read` returns. `getOneMarker` and
`posSampleInPgen` are never called on this path.

This matters because `PgenShared::readPsamFile` keys its sample map on psam
**column 0** while its comment says IID: with a `#FID IID ...` header, column 0
is FID. A new `PgenShared(pgen, psam, pvar)` overload takes the identity
mapping and **counts psam data rows directly** rather than taking
`sampleInPgen.size()` -- which would undercount whenever two samples share an
FID. See LEDGER #52 for the underlying entry; step 2 is untouched.

The R side reads IID **by name** out of the psam header
(`readGenoSampleIDs` in `R/Geno.R`), which is what makes the phenotype merge
correct regardless of whether FID == IID.

## Dosages are hard-called, once, loudly

A 2-bit store cannot hold a dosage. `PgenByteSource` rounds to the nearest of
{0,1,2}, counts how many genotypes it had to round, and raises ONE R warning at
the end of the scan naming the count. Detection is a non-integrality test on the
decoded doubles, so an integral dosage is correctly not counted.
`extdata/input/dosage_20markers_100samples.{pgen,pvar,psam}` pins this: 20
markers x 100 samples whose ALT dosage cycles 0.1 ... 1.0, giving 1,800
fractional and 200 exactly-integral genotypes.

## Gates

- `tests/testthat/test-step1-pgen.R` -- fast-tier unit tests over the R seams,
  a standard-tier comparison of the whole loaded GRM genotype store
  (`M`, `N`, `alleleFreqVec`, `MACVec`, `getQCdMarkerIndex`, 50 decoded markers,
  all `expect_identical`), the standard-tier step-0 golden, the dosage warning,
  and the full-tier five-configuration step-1 byte comparison.
- `harness/make_goldens.sh` runs step 0 and all five step-1 configs through the
  pgen reader and `cmp`s each against its PLINK twin **inside the candidate
  directory** — byte-exact, where `compare_goldens.sh`'s step-1 arm is only
  toleranced. The pgen outputs are also committed as goldens.
  **Not committing them was tried first, and the full tier caught it**:
  `compare_goldens.sh`'s step0 loop has an "EXTRA in candidate" arm that flags
  any candidate file the goldens do not know about, and the two step-0 pgen
  outputs tripped it (exactly the 2 harness failures seen). Step1 has no such
  arm, so its five would have passed unnoticed — the worse of the two failure
  modes, and the reason the fix was to commit the files rather than to silence
  the check.

## One trap, met while building this

**`createSparseGRM` is not deterministic across processes on its own.** It draws
5,000 of the 9,588 MAF-passing markers with `sample()`, and the `set.seed(1)`
lives in `extdata/createSparseGRM.R:99`, NOT in the function. A first attempt at
the step-0 comparison called `createSparseGRM()` directly from a script with no
seed; the two `.mtx` files differed (4,093 vs 4,094 nonzeros) and it looked like
a reader bug. The control settled it: **the PLINK arm did not reproduce
itself either.** Any step-0 comparison must go through the CLI, or set the seed
itself. The genotype-store comparison in the test file exists partly for this
reason -- it is RNG-free and localises a difference to the reader or not.
