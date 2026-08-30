#!/usr/bin/env bash
# Run the test suite.
#   bash harness/run_tests.sh [fast|standard|full]
#
# fast     - unit + golden-invariant + comparator self-tests only (seconds)
# standard - fast + step-2 end-to-end runs vs goldens (~3-4 min)  [default]
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
pixi run Rscript -e 'res <- testthat::test_dir("tests/testthat", stop_on_failure = TRUE)'
