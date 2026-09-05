#!/usr/bin/env bash
# Run the test suite.
#   bash harness/run_tests.sh [fast|standard|full]
#
# fast     - unit + golden-invariant + comparator self-tests only (seconds)
# standard - fast + step-2 end-to-end runs vs goldens (~2.5 min)  [default]
# full     - standard + step-0/step-1 regeneration with toleranced compare,
#            step-1 determinism double-run, and the permuted-phenotype
#            negative control (~20+ min)
#
# SAIGE_COMPARE_MODE=byte|numeric (default byte): byte is valid only on the
# platform that generated harness/goldens (see CLAUDE.md); numeric parses the
# tables and compares numerics at a relative tolerance for cross-platform CI.

set -uo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
source "$REPO/harness/env.sh"
export SAIGE_TEST_LEVEL="${1:-standard}"
export SAIGE_COMPARE_MODE="${SAIGE_COMPARE_MODE:-byte}"
echo "== SAIGE tests: level=$SAIGE_TEST_LEVEL compare=$SAIGE_COMPARE_MODE"
cd "$REPO"

# Call the pinned env's Rscript directly rather than through `pixi run`: the
# wrapper costs 2.5-5s per invocation and adds nothing here (env.sh above
# already exports the thread pins, and the R shell script sets its own R_HOME,
# verified to load SAIGE from a stripped environment).
RSCRIPT="$REPO/.pixi/envs/default/lib/R/bin/Rscript"

# Test files run in parallel (Config/testthat/{edition,parallel} in DESCRIPTION).
# Workers block on their step-2 child processes rather than spinning, so one
# worker per core is right; env.sh pins every child to a single thread, so the
# workers cannot oversubscribe each other.
export TESTTHAT_CPUS="${TESTTHAT_CPUS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"

# The ledger must parse before anything else runs: a bare `|` in a row's prose
# once hid a row from `ledger.py open` for a day (LEDGER #66), and the open
# count quoted in CLAUDE.md drifted from the tool's. Every tier, and CI, gates it.
python3 "$REPO/harness/tools/ledger.py" lint || exit 1

# package/load_package are REQUIRED for the parallel path: test_dir defaults to
# load_package="none", which hands each callr worker an empty package name and
# fails with "attempt to use zero-length variable name" before any test runs.
"$RSCRIPT" -e 'res <- testthat::test_dir("tests/testthat", package = "SAIGE", load_package = "installed", stop_on_failure = TRUE)' || exit 1
