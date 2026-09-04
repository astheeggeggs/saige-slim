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
