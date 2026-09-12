# Five open issues in the region path, and what each costs to ignore

Written 2026-09-12, after fixing LEDGER #129 and #130. **Nothing here is
urgent, nothing here is half-implemented, and every one of them is upstream's
behaviour too** (verified against `main` at 365deeff, not assumed). The tree is
clean and pushed. This page exists so the list can be read in one sitting
instead of reconstructed from six ledger rows.

## The one idea behind all five

SAIGE's single-variant code is correct for one variant's dosage vector. The
gene-based tests reuse four of its idioms -- **fold to the minor allele**,
**MAF**, **`2 - G`**, and **the dosage-zeroing filter** -- on quantities that
are not one variant's dosage vector: a bundle of rare variants squashed
together, a sum across a whole gene, a vector of user weights. Each of the five
issues below is one of those four idioms applied to one of those quantities.
Two more of the same kind were fixed today (#129, #130).

That is also why the owner's instruction -- "an allele count should be a MAC" --
is the general fix for the whole column rather than five separate repairs.

## The five

| # | What it is | Cost of ignoring it | Moves goldens if fixed? |
|---|---|---|---|
| **#134** | The squashed-variant bundle's variance-ratio band, sparse-GRM choice and exact-test eligibility are decided by the **weighted** sum, not by an allele count. So the *scale* of the weights in a group file decides which statistical treatment a set of alleles gets. | Real and live. Measured on the bundled fixture: **8 of 12** bundles land in the wrong band. `GENE1 lof` at max MAF 0.001 is 5 alleles, and it goes to the top band under the CADD weights (38.4) and the bottom band under REVEL (3.34) -- same five alleles. | **Yes** -- the multiWeights and FlexRV goldens and their PGEN twins. The only one of the five that changes numbers in a DEFAULT configuration. |
| **#135** | The burden emitter divides a cell's **aggregate** allele count by 2N and folds it as if it were a frequency, so the MAC it gates on *decreases* as the cell's allele count grows past N. | Currently **inert** -- on golden P65 the affected cells all land in the same band anyway, so no output is wrong today. Bites when the folded-down value crosses a band boundary or the min-group-MAC gate: an aggregate of 1995 at N = 1,000 is gated as MAC 5, and `--minGroupMAC_in_BurdenTest=10` would refuse the cell for being too common. | **No** -- needs a new golden with a cell whose aggregate count approaches 2N. |
| **#133** | A cell's `MAC` is the folded count while `MAC_case` and `MAC_control` are raw un-folded halves, so they do not sum to it. | Cosmetic but misleading: 444 against 151 + 1405 on P62/P63. Anyone summing the halves to check the MAC, or computing a case/control MAC ratio, gets the wrong allele. | **Yes**, but only the three goldens with cells above AF 0.5 (P62, P63, P65) -- below 0.5 the folded and un-folded counts agree. |
| **#132** | ACAT-V's combination weights mix the user's weights with a hardcoded `1` for each squashed bundle, and CCT normalises by their sum. | Multiplying every weight in a group file by a constant -- which should be a no-op, and is byte-exactly a no-op for every other p-value -- moves `Pvalue_ACATV` by up to **10.5 %** and `Pvalue_ACATO` by 3.8 %. Only in multi-weight mode, where the ACAT columns exist. | **Yes** -- the four multiWeights goldens and the FlexRV SKAT-O golden. Not the BURDEN-only ones (no ACAT columns). |
| **#131** (deferred half) | An annotated variant whose **alternate** allele is the common one stays in its annotated bin, is admitted and weighted as *rare* on its folded MAF, and the burden counts the **reference** allele -- the one the annotation does not describe. **DECIDED 2026-09-12: leave as is, for now.** | A cell labelled "rare LoF" can be driven by a rare allele that carries no annotation. How much that matters depends entirely on how often such sites occur in real data. | n/a while the decision stands. |

## The one measurement that would change a decision

**#131's revisit trigger: the fraction of annotated sites in a real callset
whose ALT allele is the major one.** It cannot be measured here -- the bundled
fixture has 0 of 70 markers above AF 0.5, and `--AlleleOrder=ref-first` inverts
all 70 rather than giving a realistic mixture. The 50k and AoU fixtures that
could have shown it were cleared from disk.

- If that fraction is negligible -- which it will be wherever a pipeline already
  filters on MAF or AAF before SAIGE sees the data -- then the decision costs
  nothing and option (c) (warn or refuse) would be a warning that never fires.
- If it is a few percent, the choice moves real results and deserves a fixture
  and goldens of its own.

## If you fix them, fix them in this order

1. **#134** -- the only one that is wrong in a default configuration, and the
   fix is to read a quantity already computed 13 lines above the three
   decisions. Its own sibling (the burden emitter) already does it right;
   LEDGER #118 fixed one of the pair and missed the other.
2. **#135** -- same root cause, one pass later, and the fix (`MAC = altCounts`)
   makes the printed and gating MACs one number and deletes a clamp.
3. **#133** -- small, contained, and makes the output self-consistent.
4. **#132** -- a different axis (weight commensurability, not the fold), and it
   needs a choice about what a squashed bundle's combination weight should be.

#134 and #135 share a root cause and a test shape (the same alleles must get
the same treatment regardless of the weights or the aggregation), so they are
better done together than apart.

## Unverified leads

The audit that found #134 and #135 was stopped after four of six surfaces. Its
remaining candidates are recorded here as LEADS -- **I did not verify these, and
an agent's framing has been wrong twice today even when its finding was right**:

- `--impute_method=minor` is documented as imputing minor-allele homozygotes;
  it may in fact impute zero copies of the minor allele, i.e. the
  major-allele homozygote (`R/SAIGE_Test_main.R:82`).
- `--is_output_moreDetails` documents four column names (`homN_Allele2_cases`
  and siblings) that no writer may emit (`R/SAIGE_Test_main.R:89`).
- The squashed bundle's weighted sum is emitted under the `AC_Allele2` /
  `AF_Allele2` column names, which name an allele it does not have
  (`src/Main.cpp:3675`).
- The FlexRV reconciliation control in `helper-step2-flexrv.R:176` may evaluate
  its transforms at the UN-folded ALT frequency. If so, the check that validates
  the FlexRV weights has been comparing against the wrong frequency -- worth
  checking first of these four, because it is a check rather than an output.

Re-running the audit: the script is under
`~/.claude/projects/.../workflows/scripts/audit-allele-convention-*.js`.
