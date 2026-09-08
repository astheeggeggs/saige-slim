# Retired genotype fixtures (bgen / vcf / sav)

51 files, 18 MB, moved here from `extdata/input/` on 2026-09-08 (LEDGER #96).

They are the inputs for the BGEN, VCF and SAV readers, which this branch
stripped in Phase 1. Nothing under `R/`, `src/`, `tests/`, `extdata/` or
`harness/` can open them, and nothing references them.

They are kept for the same reason `../step1` and `../step2` keep the retired
bgen/vcf GOLDENS: if a reader is ever restored, the inputs and the expected
outputs are both still here and still paired.

**Why they moved rather than stayed.** In `extdata/input/` they sat beside the
live fixtures under the same base names --
`genotype_100markers.{bed,pgen,bgen,vcf,vcf.gz,sav}` -- so choosing a fixture
by name could hand you one the package cannot read, and the failure surfaced as
`long flag "bgenFile" is invalid` rather than as anything about the format.

`readme-upstream-recipes` is upstream's note on how these files were generated
(bgzip/tabix/qctool/plink2/bgenix); it came from `extdata/input/readme`, where
every line was about a file in this directory.

The plink `.log` files documenting the generation runs stayed in
`extdata/input/`, because those runs also produced the live bed/bim/fam.
