# Determinism environment for all harness/golden runs. Source before ANY step0/1/2 run.
#
# Step 1 is NOT reproducible without these: the float parallel_reduce noise feeds
# back into .Random.seed via GetTrace, changing which variance-ratio markers get
# sampled.
#
# The original note here added "--nThreads=1 does not serialize TBB
# (setThreadOptions is only called when nThreads > 1,
# R/SAIGE_fitGLMM_fast.R:1055)". That is STALE: setThreadOptions is called
# UNCONDITIONALLY now (R/SAIGE_fitGLMM_fast.R:824, "LEDGER #7 fix"). Do not
# conclude from the old wording that these exports are redundant -- they still
# pin OpenBLAS and OpenMP, which that call does not touch.
#
# And note these are NOT sufficient at large N: LEDGER #56 records step 1 giving
# six distinct variance-ratio files over eleven identical runs at N = 50,000
# with --useSparseGRMforVarRatio=TRUE, WITH this file sourced and nThreads
# forced to 1. Sourcing env.sh buys reproducibility on the bundled N = 1,000
# fixtures; it does not buy it everywhere.
export RCPP_PARALLEL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# The step-1/step-2 eligibility gates (R/FitGates.R) ship ON: a binary fit with
# fewer than 100 cases, a categorical covariate level with fewer than 10 cases
# or controls, or a separated covariate GLM is REFUSED. The bundled fixture is
# below two of those by construction -- 98 cases, and a10 level 2 has 7 -- so
# every golden and every test in this repo would be refused. That refusal is
# correct (the fixture is a mechanical fixture, not a valid analysis), so the
# harness turns the gates off here rather than on ~150 command lines. A run
# with them off prints a WARNING line saying so; tests/testthat/test-step1-fit-gates.R
# sets SAIGE_FIT_GATES=1 per run and is the only place the shipped defaults
# are exercised.
export SAIGE_FIT_GATES=0
