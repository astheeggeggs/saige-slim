# Handoff 2026-09-15: LEDGER #151's guard, `harness/tools/optiondrift.py`

The ledger's five open rows are all decisions or an unreproduced crash
(#146b, #150, #151, #153, #154). The one piece of prescribed work among them
was #151's STATUS: the tool that stops the NEXT CLI-vs-API default
divergence. It exists now and runs before every tier.

## What it does

`harness/tools/optiondrift.py` (pure Python, no R, seconds) paren-matches
every `make_option(` block of the four `extdata/*.R` wrappers, the exported
function's formals (`fitNULLGLMM_multi` + `fitNULLGLMM` for step 1, the
multi's own winning on overlap), and the wrapper's call site to map flag to
formal (`opt$x`, `opt[["x"]]`, or a local whose multi-line assignment
mentions exactly one `opt$x`). Values are compared, not strings: `1L`/`1`,
`1e-05`/`1e-5`, a comma string against `c(...)`, and three coercions the
readers perform (`""`/NULL/`character(0)`; `"1"`/1 for the sex codes;
`"FALSE"`/FALSE for `--invNormalize`), each documented in `equivalent()`
with the reader line it rests on. An allowed divergence carries its ledger
row and BOTH values; a moved value fails the run as STALE, a divergence that
stops existing fails it too.

Controls, all drilled: `selftest` plants one divergence and three
representation-only agreements in a synthetic wrapper (caught / passed);
deleting an ALLOWED entry fails the run with `DRIFT`; changing one's value
fails it with `STALE`. `run_tests.sh fast` green with the wiring.

## What its first run found

Twelve divergences, not four. The eight new ones are on #151's STATUS and
the ALLOWED list. Two are LIVE in the arithmetic: step 1 `--relatednessCutoff`
0 vs formal 0.125 (the `getsubGRM` cutoff on a loaded sparse GRM), and step 0
`--numRandomMarkerforSparseKin` 2000 vs 1000 (the row had called that split
inert, which is true only of `getNeff`; the step-0 golden pins 2000 because
the CLI made it). The rest: step-0 `--nThreads` 16 vs 1,
`--sampleIDColinphenoFile` `IID` vs `""` (steps 1 and getNeff: an API caller
on the default indexes the phenotype file by an empty column name), and
`--outputPrefix` `~/` vs `""` (steps 0, 1, getNeff).

Not compared: the two FlexRV spec defaults (a named constant on both sides,
reported as `skip` under `--all`), and the `[default=...]` prose in help
strings and `@param` lines -- the `--X_PARregion` third value stays a
recorded, unguarded fact on the row.

## Gate

No `src/` or `R/` change; no tier owed. `ledger.py lint`, `docdrift.py`,
`optiondrift.py selftest`, `run_tests.sh fast` all green.

## What is still open

Every open row is the owner's: #146(b) all-zero weight matrix (recommend the
stop stands), #150 (code stays, surfaces corrected), #151 (twelve
per-divergence calls, now with a guard so the count cannot grow unnoticed),
#153 (`getNeff`'s `tau = (0, 1)`), #154 (one teardown segfault under load,
0/8 in isolation).
