# Retired genotype fixtures (bgen / vcf / sav) and their group files

51 files, 18 MB, moved here from `extdata/input/` on 2026-09-08 (LEDGER #96),
plus the SIX `group_new_chrposa1a2*.txt` group files moved on the same day
(LEDGER #113).

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


## The `group_new_chrposa1a2*.txt` group files (LEDGER #113)

`group_new_chrposa1a2{,_2chr,_emptychunk,_test,_withWeights,_withchr}.txt`,
24 KB, moved 2026-09-08. They key markers as `chr:pos:A1:A2`, which is what the
BGEN and VCF readers used; **PLINK and PGEN key on bim column 2 (rs IDs), so
every one of them matches ZERO markers in every live fixture.**

They were not merely unused, they were a trap. A group file that matches nothing
does not fail: the region run writes a ZERO-BYTE results file, exits 0 and
prints "Analysis done!" (LEDGER #113 -- the same symptom #76 was raised for, by
a route #76's `--chrom` guard does not cover). Reproduced on the bundled
fixture: `genotype_100markers` + `group_new_chrposa1a2.txt` + `--chrom=1` logs
"101 markers in 'RegionFile' are not in 'GenoFile'." and produces nothing.

Usage audit before the move: these six were referenced 0, 1 or 2 times and every
reference was a COMMENT (`extdata/cmd.sh`, `harness/make_goldens.sh`) or a docs
note that already said they were "unusable on the PLINK path"
(`docs/design/phase5-multiweights-design.md`). The live `group_new_snpid*.txt`
three are referenced 25, 9 and 7 times by goldens and tests and stayed.

Moving them does not fix #113's guard -- a hand-written group file with a typo
reaches the same silent path -- it removes the easiest way to trip it.