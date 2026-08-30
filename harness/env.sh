# Determinism environment for all harness/golden runs. Source before ANY step0/1/2 run.
#
# Step 1 is NOT reproducible without these: --nThreads=1 does not serialize TBB
# (setThreadOptions is only called when nThreads > 1, R/SAIGE_fitGLMM_fast.R:1055),
# and the float parallel_reduce noise feeds back into .Random.seed via GetTrace,
# changing which variance-ratio markers get sampled.
export RCPP_PARALLEL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
