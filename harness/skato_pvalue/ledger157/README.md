# LEDGER #157 local reproduction (2026-09-19): does not reproduce on sim50k

#157 observed, on All of Us v8 exome chr20 (N = 70,272, 9,515 cells), that a
PLINK-1 `.bed` subset to a subset of its OWN markers perturbs 4 cells at
~2.2e-07, while a PGEN written from the same `plink2` invocation reproduces the
full-bed output exactly. Mechanism not identified.

`run3.sh` is the same three-way design on sim50k (N = 50,000, 210 genes, 2,447
tested cells, `--is_fastTest=FALSE`, sparse GRM so #106's PCG amplifier is
live). It builds a bed subset and a pgen of the SAME 37,765 markers -- the exact
set named in the group file, so the subsetting is inert for the same reason it
was on AoU -- and runs the region path over all three.

**Result: all three byte-identical**, on the results file, `.singleAssoc.txt`
and `.skatoMethod.txt`, with 0 numeric differences in every pairing.

The control that makes that negative meaningful: the arms were verified to
differ. `fullbed` opened `sim50k/chr20_sim.bed` (210,515 markers), `subbed`
opened `sub_bed.bed` (37,765), `subpgen` took the PGEN reader; wall times 104 /
97 / 80 s; all three tested the identical 2,447 cells. A negative from three
arms that turned out to be the same run would prove nothing.

So the asymmetry is not a generic property of the PLINK reader, and the open
question narrows to "requires a cell on #106's iteration boundary, which this
fixture has none of". The row's cheapest next step is unchanged and still needs
All of Us.
