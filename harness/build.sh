#!/usr/bin/env bash
# Canonical build for gates: ALWAYS from clean objects.
#   bash harness/build.sh
# R's make has no header dependency tracking; a stale dependent .o after a
# class-layout change in a header produces intermittent memory corruption
# (Phase 4: two members added to PlinkClass, GenoReaderRegistry.o stale ->
# segfaults inside readFamFile). Cleaning costs ~4 min; debugging cost hours.
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO"
rm -f src/*.o src/*.so
pixi run Rscript -e 'Rcpp::compileAttributes(".")' > /dev/null
pixi run R CMD INSTALL . > harness/work/build.log 2>&1 || { echo "INSTALL FAILED"; grep -iE "error" harness/work/build.log | head -20; exit 1; }
echo "clean install OK ($(git rev-parse --short HEAD))"
