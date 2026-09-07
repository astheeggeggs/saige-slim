#!/usr/bin/env bash
# Guard Malloc over a step-1 LOCO fit -- the characterization of LEDGER #102.
#
#   bash harness/guardmalloc.sh [clean|plant|both]     (default: both)
#
# clean  run the two-chromosome LOCO fit under /usr/lib/libgmalloc.dylib. Every
#        heap allocation gets its own page with an unmapped page after it, so
#        ANY read past the end of a malloc'd buffer faults on the instruction
#        that does it. Expect exit 0.
# plant  the same fit, with one deliberate read of geno.alleleFreqVec past its
#        end injected at the same program point (Get_OneSNP_StdGeno's first act
#        is alleleFreqVec[SNPIdx], and the fixture has 2,000 markers of which
#        1,666 pass --minMAFforGRM=0.01, so index 1999 is 334 elements out).
#        Expect a fatal signal UNDER the guard and a clean exit WITHOUT it.
#
# The plant is not colour. Guard Malloc fails silently in exactly the way that
# makes a clean run look like a pass:
#
# 1. `Rscript` LOSES DYLD_INSERT_LIBRARIES. It routes through a SIP-protected
#    binary, which strips DYLD_*; the GuardMalloc banner still prints (from the
#    Rscript process, before the exec) while the R process that runs the fit
#    allocates from the ordinary system malloc. Measured 2026-09-07: with the
#    banner in the log, an arma::fvec sized 1666 read 334 elements past its end
#    and returned zeros. So the banner proves nothing; only the plant does.
#    Hence $R_HOME/bin/exec/R directly, as in harness/tsan.sh.
# 2. Bypassing bin/R also loses the DYLD_FALLBACK_LIBRARY_PATH that etc/ldpaths
#    sets, without which SAIGE.so cannot resolve libplink2_includes.a.
# 3. ARMA_NO_DEBUG is irrelevant here either way: the OOB access is
#    `alleleFreqVec[i]`, and armadillo's operator[] is unchecked in every build.
#    A bounds-checked (ARMA_NO_DEBUG-less) build does NOT catch this; only a
#    page-granular allocator does.
#
# SAIGE_GUARD_LIB points R_LIBS at another install (used to show this red at
# unpatched HEAD). SAIGE_GUARD_TIMEOUT caps each run, default 300s.
set -euo pipefail
MODE="${1:-both}"
case "$MODE" in clean|plant|both) ;; *) echo "usage: guardmalloc.sh [clean|plant|both]" >&2; exit 2 ;; esac
[ "$(uname -s)" = Darwin ] || { echo "guardmalloc.sh is macOS-only (libgmalloc)" >&2; exit 2; }
GMALLOC=/usr/lib/libgmalloc.dylib
[ -f "$GMALLOC" ] || { echo "no $GMALLOC on this machine" >&2; exit 2; }

REPO="$(cd "$(dirname "$0")/.." && pwd)"
source "$REPO/harness/timeout.sh"
PREFIX="$REPO/.pixi/envs/default"
WORK="${SAIGE_GUARD_WORK:-$REPO/harness/work/guardmalloc}"
LIMIT="${SAIGE_GUARD_TIMEOUT:-300}"
LIB="${SAIGE_GUARD_LIB:-$PREFIX/lib/R/library}"
rm -rf "$WORK"; mkdir -p "$WORK/fix"

# --- the fixture: the same two-chromosome subset test-step1-lowmem-loco.R uses,
# built through the same helper so the two cannot drift apart.
SAIGE_GUARD_REPO="$REPO" SAIGE_GUARD_FIX="$WORK/fix" "$PREFIX/lib/R/bin/Rscript" -e '
  repo <- Sys.getenv("SAIGE_GUARD_REPO")
  saige_extdata <- file.path(repo, "extdata")
  source(file.path(repo, "tests/testthat/helper-pheno-scale.R"))
  sm_two_chr_fixture(file.path(Sys.getenv("SAIGE_GUARD_FIX"), "sub"))' \
  > "$WORK/fixture.log" 2>&1

# The plant driver: library, then one out-of-range read at the LOCO entry point.
# Survival is signalled by CREATING A FILE, not by printing: R deparses a
# tracer into its segfault traceback, so any sentinel string in the expression
# appears in the log of the run that CRASHED at it.
cat > "$WORK/plant.R" <<'PLANT'
suppressMessages(library(SAIGE))
trace(setStartEndIndex,
      tracer = quote({
        invisible(Get_OneSNP_StdGeno(1999L))
        file.create(Sys.getenv("SAIGE_GUARD_SENTINEL"))
      }),
      where = asNamespace("SAIGE"), print = FALSE)
source("step1_fitNULLGLMM.R")
PLANT

# --- the launcher. See notes 1 and 2 in the header for why it is spelled this way.
cat > "$WORK/run_r.sh" <<EOF
#!/bin/bash
export R_HOME="$PREFIX/lib/R"
export DYLD_FALLBACK_LIBRARY_PATH="\$R_HOME/lib:$PREFIX/lib:\$HOME/lib:/usr/local/lib:/usr/lib"
export R_LIBS="$LIB:\$R_HOME/library"
export RCPP_PARALLEL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
[ -n "\${SAIGE_GUARD_ON:-}" ] && export DYLD_INSERT_LIBRARIES="$GMALLOC"
exec "\$R_HOME/bin/exec/R" --no-save --no-restore --no-echo --no-init-file "\$@"
EOF
chmod +x "$WORK/run_r.sh"

FAILED=0
gm_run() { # gm_run <name> <driver> <guard on|off>; leaves the status in $RWT_STATUS
  local name="$1" driver="$2" guard="$3"
  mkdir -p "$WORK/out.$name"; rm -f "$WORK/$name.survived"
  ( cd "$REPO/extdata" && SAIGE_GUARD_ON="$guard" SAIGE_GUARD_SENTINEL="$WORK/$name.survived" \
    run_with_timeout "$LIMIT" "$WORK/$name.log" \
      bash "$WORK/run_r.sh" --file="$driver" --args \
      --plinkFile="$WORK/fix/sub" --nThreads=1 --minMAFforGRM=0.01 --LOCO=TRUE \
      --isLowMemLOCO=TRUE --skipVarianceRatioEstimation=TRUE \
      --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
      --phenoCol=y_binary --covarColList=x1,x2,a9,a10 --qCovarColList=a9,a10 \
      --sampleIDColinphenoFile=IID --traitType=binary \
      --outputPrefix="$WORK/out.$name/fit"
    echo "$RWT_STATUS" > "$WORK/$name.rc" )
  RWT_STATUS=$(cat "$WORK/$name.rc")
}

verdict() { # verdict <label> <expected: ok|dead> <status>
  local label="$1" want="$2" st="$3" got
  case "$st" in 0) got=ok ;; 124) got="TIMED OUT" ;; *) got=dead ;; esac
  if [ "$got" = "$want" ]; then
    printf '  %-38s status=%-3s %s\n' "$label" "$st" "PASS"
  else
    printf '  %-38s status=%-3s %s (wanted %s, got %s)\n' "$label" "$st" "FAIL" "$want" "$got"
    FAILED=1
  fi
}

echo "=== Guard Malloc ($MODE); R_LIBS=$LIB"
if [ "$MODE" = clean ] || [ "$MODE" = both ]; then
  gm_run clean "$REPO/extdata/step1_fitNULLGLMM.R" on
  verdict "clean fit, guard ON" ok "$RWT_STATUS"
fi
if [ "$MODE" = plant ] || [ "$MODE" = both ]; then
  gm_run plant_off "$WORK/plant.R" ""
  verdict "planted OOB read, guard OFF" ok "$RWT_STATUS"
  if [ ! -f "$WORK/plant_off.survived" ]; then
    echo "  the plant never reached its read with the guard off -- the ON arm gates nothing"; FAILED=1
  fi
  gm_run plant_on "$WORK/plant.R" on
  verdict "planted OOB read, guard ON" dead "$RWT_STATUS"
  if [ -f "$WORK/plant_on.survived" ]; then
    echo "  the guard did not catch the planted read -- a clean pass proves NOTHING"; FAILED=1
  fi
fi
echo
echo "Expected: the clean fit exits 0 under the guard, the planted read survives"
echo "without it and is killed by it. logs: $WORK"
exit $FAILED
