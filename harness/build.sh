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
# builds are a no-op. Mirrors docker/Dockerfile's layer, but uses R's own
# configured C++ compiler (`R CMD config CXX`) so it works natively on macOS
# and Linux without assuming pixi exports $CXX -- it doesn't on a bare Linux
# runner (no gxx_linux-64/clangxx-style compiler package in pixi.toml), which
# is exactly what left `"$CXX"` expanding to an empty string and failing with
# "command not found" the first time this ran in CI. R's compiler is
# guaranteed to exist: this package can't build without it either way.
# TWO separate guards, deliberately. The compiled archive and the header tree are
# needed by different things and can exist independently: the archive is what the
# package links against, while the HEADERS are needed by every compile of this
# package, because src/Makevars carries -I../plink-ng/2.0/include. CI caches only
# the archive (plink-ng/ is gitignored), so a cache hit used to restore the .a,
# skip this whole block, and then fail the package build with
# "pgenlib_read.h: No such file or directory" -- correct on a cold cache, broken
# on every warm one. Guarding them separately also keeps the cache worth having:
# a warm run re-fetches the source (seconds) but skips the compile (minutes).
if [ ! -f plink-ng/2.0/include/pgenlib_read.h ]; then
  echo "plink-ng headers not found; fetching ..."
  rm -rf plink-ng plink-ng-2.0.0-a.6.16
  curl -L https://github.com/chrchang/plink-ng/archive/refs/tags/v2.0.0-a.6.16.tar.gz | tar -zx
  mv plink-ng-2.0.0-a.6.16 plink-ng
fi

if [ ! -f .pixi/envs/default/lib/libplink2_includes.a ]; then
  echo "libplink2_includes.a not found; building plink-ng ..."
  # R CMD config CXX returns "<compiler> -std=gnu++NN"; take just the
  # compiler binary so it doesn't collide with our own -std=c++14 below.
  RCXX="$(pixi run R CMD config CXX | awk '{print $1}')"
  pixi run bash -c "'$RCXX' -std=c++14 -fPIC -O3 -I.pixi/envs/default/include -L.pixi/envs/default/lib -o libplink2_includes.a plink-ng/2.0/include/*.cc -shared -lz -lzstd -lpthread -lm -ldeflate"
  mv libplink2_includes.a .pixi/envs/default/lib
fi

rm -f src/*.o src/*.so
pixi run Rscript -e 'Rcpp::compileAttributes(".")' > /dev/null
mkdir -p harness/work
pixi run R CMD INSTALL . > harness/work/build.log 2>&1 || { echo "INSTALL FAILED"; grep -iE "error" harness/work/build.log | head -20; exit 1; }
echo "clean install OK ($(git rev-parse --short HEAD))"
