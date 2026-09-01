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

# PGEN support (Phase 7) needs pgenlib from plink-ng, compiled into
# .pixi/envs/default/lib/libplink2_includes.a. Not vendored in git
# (.gitignore excludes plink-ng/); fetch+compile once, guarded so repeat
# builds are a no-op. Mirrors docker/Dockerfile's layer, but uses the
# pixi-env compiler ($CXX) so it works natively on macOS and Linux, not just
# the Linux Docker image's fixed cross-compiler name.
if [ ! -f .pixi/envs/default/lib/libplink2_includes.a ]; then
  echo "libplink2_includes.a not found; fetching + building plink-ng ..."
  if [ ! -d plink-ng ]; then
    curl -L https://github.com/chrchang/plink-ng/archive/refs/tags/v2.0.0-a.6.16.tar.gz | tar -zx
    mv plink-ng-2.0.0-a.6.16 plink-ng
  fi
  pixi run bash -c '"$CXX" -std=c++14 -fPIC -O3 -I.pixi/envs/default/include -L.pixi/envs/default/lib -o libplink2_includes.a plink-ng/2.0/include/*.cc -shared -lz -lzstd -lpthread -lm -ldeflate'
  mv libplink2_includes.a .pixi/envs/default/lib
fi

rm -f src/*.o src/*.so
pixi run Rscript -e 'Rcpp::compileAttributes(".")' > /dev/null
pixi run R CMD INSTALL . > harness/work/build.log 2>&1 || { echo "INSTALL FAILED"; grep -iE "error" harness/work/build.log | head -20; exit 1; }
echo "clean install OK ($(git rev-parse --short HEAD))"
