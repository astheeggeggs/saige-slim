#!/usr/bin/env bash

# osx-arm64 FLIBS TRAP (2026-09-19, first build ever attempted on Apple Silicon).
# `R CMD INSTALL` fails at the link step with:
#     ld: warning: directory not found for option
#         '-L<env>/lib/gcc/arm64-apple-darwin20.0.0/14.3.0'
#     ld: library not found for -lheapt_w
# CAUSE: the osx-arm64 r-base in pixi.lock baked FLIBS against
# libgfortran-devel 14.3.0, while the lock installs 14.4.0. Only the VERSION
# component of the -L path is wrong -- libheapt_w.a is present, at
# lib/gcc/arm64-apple-darwin20.0.0/14.4.0/. So this is a -L fix, not a missing
# library, and NOT the same symptom as the x86_64 case (where the baked
# directory happens to exist).
# WORKAROUND, local and non-invasive -- do NOT edit pixi.toml or pixi.lock,
# which are deliberately pinned and whose regeneration is its own gated change:
#     printf 'FLIBS = -L%s/lib/gcc/arm64-apple-darwin20.0.0/14.4.0 -lheapt_w -lgfortran -lquadmath -lm\n' \
#       "$PWD/.pixi/envs/default" > /tmp/saige-makevars
#     R_MAKEVARS_USER=/tmp/saige-makevars bash harness/build.sh
# (R_MAKEVARS_USER rather than ~/.R/Makevars: this is per-env, not per-user.)
# Verified: `clean install OK` with the warnings ratchet still at baseline, so
# only the link line moved. osx-arm64 is declared in pixi.toml's platforms but
# no gate has ever built it -- CI is linux-64 and local dev was x86_64.

# Canonical build for gates: ALWAYS from clean objects.
#   bash harness/build.sh                 clean (the default, and what a gate wants)
#   bash harness/build.sh --incremental   keep src/*.o; REFUSES if a header moved
#
# R's make has no header dependency tracking; a stale dependent .o after a
# class-layout change in a header produces intermittent memory corruption
# (Phase 4: two members added to PlinkClass, GenoReaderRegistry.o stale ->
# segfaults inside readFamFile). Cleaning costs ~4 min; debugging cost hours.
#
# --incremental exists because most edits are to ONE .cpp and a clean rebuild of
# all 20 objects is then ~4 minutes of waiting for nothing. It does not ASK you
# to remember the clean-objects rule, it ENFORCES it: no .hpp or .h may be newer
# than the OLDEST object, and the oldest is the right comparison because a
# single object compiled before the header changed is enough to corrupt the
# build. If any is newer it refuses and names both files. With no objects
# present at all, incremental IS clean, so it proceeds.
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

INCREMENTAL=0
for arg in "$@"; do
  case "$arg" in
    --incremental) INCREMENTAL=1 ;;
    *) echo "build.sh: unknown argument '$arg' (only --incremental)" >&2; exit 2 ;;
  esac
done

if [ "$INCREMENTAL" = 1 ]; then
  # The clean-objects rule, enforced. `ls -tr | head -1` is the OLDEST object:
  # comparing against the newest would pass a tree where one straggler was
  # compiled before the header change, which is exactly the failure mode.
  oldest_o=$(ls -tr src/*.o 2>/dev/null | head -1 || true)
  if [ -n "$oldest_o" ]; then
    newer_hdr=$(find src \( -name '*.hpp' -o -name '*.h' \) -newer "$oldest_o" | head -5)
    if [ -n "$newer_hdr" ]; then
      echo "build.sh: --incremental REFUSED -- a header is newer than the oldest object." >&2
      echo "  oldest object:  $oldest_o" >&2
      echo "$newer_hdr" | sed 's/^/  newer header:  /' >&2
      echo "  R's make has no header dependency tracking, so every object compiled" >&2
      echo "  before that header is silently stale, and a stale dependent after a" >&2
      echo "  class-layout change gives intermittent memory corruption." >&2
      echo "  Run: bash harness/build.sh   (clean)" >&2
      exit 1
    fi
  fi
else
  rm -f src/*.o src/*.so
fi

pixi run Rscript -e 'Rcpp::compileAttributes(".")' > /dev/null
mkdir -p harness/work
# On failure print the DIAGNOSTICS, not every line containing "error": the
# old `grep -iE "error" | head -20` was filled by Rcpp's `error_occured`
# symbol and CI run 34837530220 failed with the real compiler error never
# shown (the log is not uploaded from a runner). "error:" with the colon is
# the compiler's and linker's spelling; "Error" R's.
# SAIGE_R_LIBRARY: install into that library instead of the default env's
# (docker/Dockerfile installs into the RUNTIME environment's R library so the
# runtime stage can take the tree whole). Empty = the default, unchanged.
pixi run R CMD INSTALL ${SAIGE_R_LIBRARY:+--library="$SAIGE_R_LIBRARY"} . > harness/work/build.log 2>&1 || {
  echo "INSTALL FAILED -- diagnostics from harness/work/build.log:"
  grep -nE "error:|^Error|ERROR:|undefined reference|Undefined symbols" harness/work/build.log | head -40
  echo "--- last 15 lines:"; tail -15 harness/work/build.log
  exit 1
}
if [ "$INCREMENTAL" = 1 ]; then
  echo "incremental install OK ($(git rev-parse --short HEAD)) -- objects kept, no header moved"
else
  echo "clean install OK ($(git rev-parse --short HEAD))"
fi

# --- the warning ratchet (2026-09-14, audit Phase D) --------------------------
# src/Makevars carried -w, which cancelled every diagnostic including the
# -Wall R passes by default, and nothing read a compiler warning anywhere. The
# goldens are a byte oracle over one fixture and cannot see an uninitialised
# read on a path the fixture never takes; the compiler can. So: count the
# src-anchored warnings (generated RcppExports.cpp excluded -- its
# cast-function-type-mismatch noise is Rcpp's, not ours) and refuse a build
# that has MORE than harness/warnings-baseline.txt records. Fewer is progress
# and is reported so the baseline can be lowered. Clean builds only: an
# incremental build compiles a subset and its count means nothing. Local clang
# only: the spellings and counts differ under CI's GCC, which is why CI sets
# SAIGE_WARNINGS_RATCHET=0 and promotes error CLASSES instead (see
# .github/workflows/tests.yml and src/Makevars).
BASELINE="harness/warnings-baseline.txt"
warn_hist() {  # src-anchored warnings from the build log, "N -Wcategory" lines
  grep 'warning:' harness/work/build.log     | grep -E '^[A-Za-z_0-9]+\.(cpp|hpp|h):' | grep -v '^RcppExports\.cpp'     | sed -E 's/.*\[(-W[^]]+)\].*/\1/' | sort | uniq -c | sort -rn | sed 's/^ *//' || true
}
if [ "${SAIGE_WARNINGS_RATCHET:-1}" = 1 ] && [ "$INCREMENTAL" != 1 ]; then
  count=$(warn_hist | awk '{s += $1} END {print s + 0}')
  if [ ! -f "$BASELINE" ]; then
    echo "build.sh: no $BASELINE; writing one at $count warnings (commit it)"
    { echo "$count"; echo "# src-anchored compiler warnings at the last accepted build, RcppExports.cpp"; echo "# excluded; first line is the ratchet. Histogram for humans:"; warn_hist | sed 's/^/# /'; } > "$BASELINE"
  else
    base=$(head -1 "$BASELINE")
    if [ "$count" -gt "$base" ]; then
      echo "build.sh: WARNING RATCHET FAILED -- $count src-anchored warnings, baseline $base." >&2
      echo "  New diagnostics (histogram now vs baseline):" >&2
      warn_hist | sed 's/^/    now:  /' >&2
      grep '^# [0-9]' "$BASELINE" | sed 's/^# /    was:  /' >&2
      echo "  Fix the new warning, or if it is deliberate raise $BASELINE in the same commit and say why." >&2
      exit 1
    elif [ "$count" -lt "$base" ]; then
      echo "warnings: $count src-anchored (baseline $base) -- fewer; lower the first line of $BASELINE to $count"
    else
      echo "warnings: $count src-anchored, at the baseline"
    fi
  fi
fi
