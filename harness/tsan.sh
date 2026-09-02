#!/usr/bin/env bash
# ThreadSanitizer over the step-2 in-process marker scan (S6).
#
#   bash harness/tsan.sh [clean|plant]     (default: clean)
#
# clean  build HEAD's sources and run the three marker workloads at 1/4/8
#        threads. Expect ZERO races and byte-identity with the serial arm.
# plant  build the SAME tree with one defect planted -- every worker shares one
#        MarkerScratch instead of owning one -- and expect races. This is not
#        optional colour: a clean TSAN run means nothing on its own, because a
#        TSAN build that sees none of this code also reports zero. Run both.
#
# Four things about this are macOS-specific and each of them fails in a way that
# looks like "no races found", which is why they are scripted rather than
# described:
#
# 1. TSAN MUST BE PRELOADED, not merely linked. The .so links the runtime, but R
#    dlopens it, by which time the interceptors cannot be installed -- the
#    runtime says so ("Interceptors are not working") and aborts.
# 2. DYLD_INSERT_LIBRARIES CANNOT BE SET FROM OUTSIDE. SIP strips DYLD_* when
#    exec'ing a protected binary, and both `Rscript` and `R` are shell wrappers
#    that route through /bin/sh. So it is exported inside an already-running
#    shell which then execs $R_HOME/bin/exec/R directly -- an unprotected binary
#    in the pixi env, which inherits it.
# 3. Bypassing R's bin/R wrapper also loses the DYLD_FALLBACK_LIBRARY_PATH that
#    etc/ldpaths sets, without which SAIGE.so's load command for
#    libplink2_includes.a resolves against nothing. Set here instead.
# 4. `R CMD INSTALL` test-loads the .so as its last step, with no preload, so it
#    aborts and DELETES the install it just made. Hence --no-test-load.
#
# The build goes to a COPY of the tree, never to src/: an in-place TSAN build
# would leave -fsanitize=thread objects behind for the next R CMD INSTALL to
# link against (the same class of hazard as the header rule in CLAUDE.md).
set -euo pipefail
MODE="${1:-clean}"
case "$MODE" in clean|plant) ;; *) echo "usage: tsan.sh [clean|plant]" >&2; exit 2 ;; esac
REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${SAIGE_TSAN_WORK:-$REPO/harness/work/tsan}"
PKG="$WORK/pkg"; LIB="$WORK/lib"; RUNS="$WORK/runs.$MODE"
PREFIX="$REPO/.pixi/envs/default"
TSAN_DYLIB="$PREFIX/lib/libclang_rt.tsan_osx_dynamic.dylib"
[ -f "$TSAN_DYLIB" ] || { echo "no TSAN runtime at $TSAN_DYLIB" >&2; exit 1; }
rm -rf "$RUNS"; mkdir -p "$WORK" "$LIB" "$RUNS"

# --- the package copy -------------------------------------------------------
rm -rf "$PKG"; mkdir -p "$PKG"
# bsdtar wants the mode flag before -C, so spell every flag out rather than
# bundling: `tar -C dir cf -` is rejected outright ("Must specify one of ...").
tar -c -f - -C "$REPO" --exclude=.git --exclude=.pixi --exclude=plink-ng \
    --exclude=harness/work --exclude='src/*.o' --exclude='src/*.so' . \
  | tar -x -f - -C "$PKG"
ln -sf "$REPO/.pixi" "$PKG/.pixi"
ln -sf "$REPO/plink-ng" "$PKG/plink-ng"

if [ "$MODE" = plant ]; then
  # The S6 drill, as a compile-time defect: hand every worker the SAME scratch.
  # It must produce worker-vs-worker races AND a wrong answer; if it produces
  # only the second, TSAN is not seeing this code and a clean run proves nothing.
  python3 - "$PKG/src/Main.cpp" <<'PY'
import sys
p = sys.argv[1]; s = open(p).read()
a = """    for(int i = 0; i < t_nStates; i++){
      m_free.push_back(make());
    }
  }"""
b = """    for(int i = 0; i < t_nStates; i++){
      m_free.push_back(make());
    }
    m_plantShared = &m_free.front()->scratch;
  }
  Marker::MarkerScratch* m_plantShared = nullptr;"""
c = "state->scratch, m_resultVec.at(i));"
d = "*m_pool.m_plantShared, m_resultVec.at(i));"
assert s.count(a) == 1 and s.count(c) == 1, "plant anchors drifted; fix harness/tsan.sh"
open(p, 'w').write(s.replace(a, b).replace(c, d))
PY
fi

# --- the instrumented build -------------------------------------------------
# -O1 -g: TSAN needs frame pointers and symbols to be readable, and -O2 inlines
# the loop body into something the reports cannot be attributed with.
# SAIGE_TSAN_ANNOTATIONS turns on the parallelFor join annotation in Main.cpp;
# see the comment there for why an uninstrumented libtbb makes it necessary.
cat > "$WORK/Makevars.tsan" <<EOF
CXX14FLAGS = -march=core2 -mtune=haswell -mssse3 -fPIC -O1 -g -fno-omit-frame-pointer -pipe -stdlib=libc++ -fvisibility-inlines-hidden -fmessage-length=0 -isystem $PREFIX/include -fsanitize=thread -DSAIGE_TSAN_ANNOTATIONS
SHLIB_CXX14LDFLAGS = -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -fsanitize=thread
EOF
echo "=== building instrumented SAIGE ($MODE) -> $LIB"
R_MAKEVARS_USER="$WORK/Makevars.tsan" pixi run --manifest-path "$REPO/pixi.toml" \
  R CMD INSTALL --no-test-load -l "$LIB" --no-docs --no-byte-compile "$PKG" > "$WORK/build.$MODE.log" 2>&1 || {
    echo "build failed; see $WORK/build.$MODE.log"; tail -20 "$WORK/build.$MODE.log"; exit 1; }

cat > "$WORK/run_r.sh" <<EOF
#!/bin/bash
export R_HOME="$PREFIX/lib/R"
export DYLD_INSERT_LIBRARIES="$TSAN_DYLIB"
export DYLD_FALLBACK_LIBRARY_PATH="\$R_HOME/lib:$PREFIX/lib:\$HOME/lib:/usr/local/lib:/usr/lib"
export R_LIBS="$LIB:\$R_HOME/library"
exec "\$R_HOME/bin/exec/R" --no-save --no-restore --no-echo --no-init-file "\$@"
EOF
chmod +x "$WORK/run_r.sh"

# --- the runs ---------------------------------------------------------------
cd "$REPO/extdata"
GOLD="$REPO/harness/goldens"
P=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr
SG="--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"

# Three workloads, for the reason the retirement criterion gives: the binary
# ones defer every MAC <= max_MAC_for_ER marker to the main thread (Rf_lchoose
# is main-thread-only), so a binary run threads LESS of itself than a
# quantitative one, and quant is where the widest concurrency actually happens.
tsan_run() { # tsan_run <workload> <threads>
  local wl="$1" nt="$2" model args
  case "$wl" in
    binary)        model=binary_sparseGRM_a;  args="" ;;
    binary_sparse) model=binary_sparseGRM_vr; args="$SG --is_fastTest=FALSE" ;;
    quant)         model=quant_sparseGRM;     args="" ;;
  esac
  local base="$RUNS/$wl.t$nt"
  SAIGE_MARKER_THREADS="$nt" \
  TSAN_OPTIONS="halt_on_error=0 history_size=7 exitcode=0 log_path=$base.tsan" \
  RCPP_PARALLEL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    pixi run --manifest-path "$REPO/pixi.toml" bash "$WORK/run_r.sh" -f step2_SPAtests.R --args \
    --bedFile=$P.bed --bimFile=$P.bim --famFile=$P.fam --AlleleOrder=alt-first \
    --minMAF=0 --minMAC=1 --LOCO=FALSE --markers_per_chunk=10000 $args \
    --GMMATmodelFile="$GOLD/step1/$model.rda" \
    --varianceRatioFile="$GOLD/step1/$model.varianceRatio.txt" \
    --SAIGEOutputFile="$base.txt" > "$base.log" 2>&1 || true
  local races threaded
  races=$(cat "$base".tsan.[0-9]* 2>/dev/null | grep -c "WARNING: ThreadSanitizer: data race" || true)
  threaded=$(grep -c "in-process threads" "$base.log" || true)
  # A threaded arm that quietly ran serially reports zero races and proves
  # nothing -- the LEDGER #30 shape, and the reason the standing test asserts
  # this same line.
  if [ "$nt" -gt 1 ] && [ "$threaded" -eq 0 ]; then
    echo "  $wl t=$nt: DID NOT THREAD (no thread-count line) -- result is meaningless"; return
  fi
  local ident="n/a"
  if [ "$nt" -gt 1 ]; then
    if cmp -s "$RUNS/$wl.t1.txt" "$base.txt"; then ident="identical"; else ident="DIFFERS"; fi
  fi
  printf '  %-14s t=%-2s races=%-4s bytes-vs-serial=%s\n' "$wl" "$nt" "$races" "$ident"
}

echo "=== TSAN runs ($MODE); 128,868 markers each"
for wl in binary binary_sparse quant; do
  for nt in 1 4 8; do tsan_run "$wl" "$nt"; done
done
echo
if [ "$MODE" = clean ]; then
  echo "Expected: races=0 everywhere and bytes-vs-serial=identical."
else
  echo "Expected: races>0 at t>1 (worker-vs-worker) and bytes-vs-serial=DIFFERS."
  echo "If this mode reports races=0, the CLEAN result above gates nothing."
fi
echo "reports and outputs: $RUNS"
