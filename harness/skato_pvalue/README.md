# SKAT-O p-value methods: the 2026-09-17 measurement harness

Results and decision: `docs/state/handoff-2026-09-17-skato-pvalue-methods.md`.
Everything here runs against the sim50k fixture (`SAIGE-sim/harness/work/sim50k`,
not in git) and the pgen reader, 1 thread, `source harness/env.sh` first.

- `profile_region.sh <scratch> skato|burden` -- the `SAIGE_PROFILE_REGION=1` run
  that produced the C++ bucket table (SKAT-O vs BURDEN, same recipe).
- `time_liu_vs_davies.R <out> <records.tsv> <sim50k dir>` -- patches
  `SKAT:::Met_SKAT_Get_Pvalue` in a child so every cell is timed under
  `optimal.adj` (Davies) and `liu.mod`; the run keeps the Davies answer and its
  output must stay byte-identical to an unpatched run (the control).
- `dump_cells.R <out> <cells.rds> <sim50k dir>` -- same patch, but records
  every cell's (Score, Phi, r.corr, SKAT p). Edit the model name for the trait.
- `variants.R <scratch> <celldir>` -- offline: reimplements SKAT's SKAT-O path
  with a pluggable inner survival function; `repro` must equal SKAT's p exactly
  (the control on the reimplementation); then modified Liu, saddlepoint
  (Kuonen), Davies acc 1e-9/lim 1e5, the split-integral ("direct") form, and
  the adaptive-Davies reference. `report.R <scratch> <celldir>` tabulates.
The Ruben positive-series reference in `variants.R` is exact but only converges
on well-conditioned spectra; on real cells it never did (SKAT keeps eigenvalues
down to mean/1e5), so it is validated on synthetic spectra only.
