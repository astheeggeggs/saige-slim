#!/usr/bin/env bash
# ThreadSanitizer over the step-2 in-process marker scan (S6).
#
#   bash harness/tsan.sh [clean|plant|plantreader]     (default: clean)
#
# clean       build HEAD's sources and run the marker workloads at 1/4/8
#             threads. Expect ZERO races and byte-identity with the serial arm.
# plant       build the SAME tree with one defect planted -- every worker shares
#             one MarkerScratch instead of owning one -- and expect races. This
#             is not optional colour: a clean TSAN run means nothing on its own,
#             because a TSAN build that sees none of this code also reports
#             zero. Run both.
# plantreader the control for the PGEN clone (S8): PgenReader::clone() hands
#             back a handle onto the SAME PgenClass instead of building a second
#             pgenlib reader, which is exactly the arrangement the split of
#             PgenShared/PgenClass exists to avoid. Run it with
#             SAIGE_TSAN_WORKLOADS=pgen; on a PLINK workload it plants nothing
#             (PLINK clones share one reader BY DESIGN -- its reads are const).
#
# SAIGE_TSAN_TIMEOUT caps each ARM (default 3600s). It exists because a plant is
# licensed to misbehave arbitrarily: on 2026-09-07 the planted `region t=4` arm
# segfaulted, R's segfault handler printed its banner and then blocked on the
# interactive abort prompt, and this script sat on it for 4h32m before it was
# killed by hand. An arm that dies or times out is reported as its own outcome,
# distinct from races=0 -- under `plant` that is the mode working, under `clean`
# it is a failure.
#
# SAIGE_TSAN_MAX_MARKERS caps every arm at the first N markers of the fixture,
# via an --idstoIncludeFile written here. It defaults to 2000 in `plantreader`
# and to the whole fixture otherwise, because the planted defect is not merely
# incorrect, it is PATHOLOGICALLY SLOW: several threads sharing one pgenlib
# reader thrash its LD cache, and the full 128,868-marker arm did not finish
# chunk 1 of 22 in eighteen minutes (the clean build does all 22 in ~2). A
# control nobody can run to completion is not a control, and the cap applies to
# the 1-thread arm too, so the byte comparison between arms stays honest.
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
case "$MODE" in clean|plant|plantreader) ;; *) echo "usage: tsan.sh [clean|plant|plantreader]" >&2; exit 2 ;; esac
REPO="$(cd "$(dirname "$0")/.." && pwd)"
source "$REPO/harness/timeout.sh"
LIMIT="${SAIGE_TSAN_TIMEOUT:-3600}"
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
  # Phase G2 (2026-09-14): the two region Workers live in src/RegionLoop.hpp
  # now, so anchors e/g are patched there and a/c in Main.cpp.
  python3 - "$PKG/src/Main.cpp" "$PKG/src/RegionLoop.hpp" <<'PY'
import sys
p = sys.argv[1]; s = open(p).read()
q = sys.argv[2]; t = open(q).read()
a = """    for(int i = 0; i < t_nStates; i++){
      m_free.push_back(make());
    }
  }"""
b = """    for(int i = 0; i < t_nStates; i++){
      m_free.push_back(make());
    }
    m_plantShared = m_free.front()->scratch[0].get();
  }
  Marker::MarkerScratch* m_plantShared = nullptr;"""
# M1: a worker state now holds one scratch PER SAMPLE-SET GROUP, so the defect
# is "every worker uses the first state's group-0 scratch". On a single-group
# workload -- which every workload here is -- that is exactly the old plant.
c = "*state->scratch[g], rows);"
d = "*m_pool.m_plantShared, rows);"
# The region loop's worker (2026-09-06) takes the same defect: every range
# prepares and tests its markers in the first state's scratch.
e = "m_fn(m_iStart + k, k, *state->scratch[0], *state->readers[0]);"
f = "m_fn(m_iStart + k, k, *m_pool.m_plantShared, *state->readers[0]);"
# The per-cell test worker of passes B and D (step 3b, 2026-09-06): same
# defect, every cell's test in the first state's scratch.
g = "m_fn((unsigned int)k, *state->scratch[0]);"
h = "m_fn((unsigned int)k, *m_pool.m_plantShared);"
assert s.count(a) == 1 and s.count(c) == 1 and t.count(e) == 1 and t.count(g) == 1, "plant anchors drifted; fix harness/tsan.sh"
open(p, 'w').write(s.replace(a, b).replace(c, d))
open(q, 'w').write(t.replace(e, f).replace(g, h))
PY
fi

if [ "$MODE" = plantreader ]; then
  # Share one pgenlib decoder across every worker: same FILE*, same fread
  # buffer, same LD cache. If S8's split is doing anything, this must produce
  # worker-vs-worker races INSIDE pgenlib and a wrong answer.
  python3 - "$PKG/src/GenoReaderRegistry.cpp" <<'PZ'
import sys
p = sys.argv[1]; s = open(p).read()
a = """    return std::unique_ptr<GenoReader>(new PgenReader(
        std::make_shared<PGEN::PgenClass>(impl().shared())));"""
b = """    return std::unique_ptr<GenoReader>(new PgenReader(m_impl));"""
assert s.count(a) == 1, "plantreader anchor drifted; fix harness/tsan.sh"
open(p, 'w').write(s.replace(a, b))
PZ
fi

# --- instrumented pgenlib ---------------------------------------------------
# pgenlib arrives as libplink2_includes.a, which harness/build.sh compiles
# WITHOUT instrumentation -- so TSAN cannot see inside the decoder, and two
# workers racing in there would be reported by nothing. That is the same blind
# spot uninstrumented libtbb has, and it matters more here: PGEN's clone (S8)
# rests entirely on the claim that separate plink2::PgenReaders over ONE
# PgenFileInfo share nothing mutable, and a clean run over an uninstrumented
# decoder does not test that claim at all. Built once and cached in $WORK.
PGENLIB_TSAN="$WORK/libplink2_includes.a"
if [ ! -f "$PGENLIB_TSAN" ]; then
  echo "=== building instrumented pgenlib -> $PGENLIB_TSAN"
  RCXX="$(pixi run --manifest-path "$REPO/pixi.toml" R CMD config CXX | awk '{print $1}')"
  (cd "$REPO" && pixi run --manifest-path "$REPO/pixi.toml" bash -c "'$RCXX' -std=c++14 -fPIC -O1 -g -fno-omit-frame-pointer -fsanitize=thread -I.pixi/envs/default/include -L.pixi/envs/default/lib -o '$PGENLIB_TSAN' plink-ng/2.0/include/*.cc -shared -lz -lzstd -lpthread -lm -ldeflate -Wl,-install_name,'$PGENLIB_TSAN'") \
    > "$WORK/build.pgenlib.log" 2>&1 || {
      echo "instrumented pgenlib build failed; see $WORK/build.pgenlib.log"; exit 1; }
fi
# NAME THE FILE, do not go through -l, and the reason is a silent failure worth
# knowing about: R puts `-L<pixi env>/lib` on the link line BEFORE PKG_LIBS, so
# `-lplink2_includes` finds the uninstrumented archive there no matter where an
# extra -L is prepended. The build then succeeds, the run works, and TSAN simply
# sees nothing inside pgenlib. Linking the absolute path in the package COPY (the
# same place the plant modes patch sources) is unambiguous, and the absolute
# install_name above makes dyld load that same file at runtime.
# (Check with `otool -L .../SAIGE.so | grep plink2`: it must print the absolute
# path under harness/work/tsan, not the bare filename.)
python3 - "$PKG/src/Makevars" "$PGENLIB_TSAN" <<'PM'
import os, sys
p, lib = sys.argv[1], sys.argv[2]
s = open(p).read()
a = "PKG_LIBS += -L.. -lplink2_includes"
b = "PKG_LIBS += %s" % lib
assert s.count(a) == 1, "Makevars anchor drifted; fix harness/tsan.sh"
open(p, 'w').write(s.replace(a, b))
PM

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
# Marker cap (see the header). The include file is written from the .pvar/.bim of
# whichever fixture a workload names, so it is built inside tsan_run.
MAXMK="${SAIGE_TSAN_MAX_MARKERS:-}"
if [ -z "$MAXMK" ] && [ "$MODE" = plantreader ]; then MAXMK=2000; fi
cd "$REPO/extdata"
GOLD="$REPO/harness/goldens"
P=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr
SG="--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"

# Three workloads, for the reason the retirement criterion gives: the binary
# ones defer every MAC <= max_MAC_for_ER marker to the main thread (Rf_lchoose
# is main-thread-only), so a binary run threads LESS of itself than a
# quantitative one, and quant is where the widest concurrency actually happens.
tsan_run() { # tsan_run <workload> <threads>
  local wl="$1" nt="$2" model args bed="--bedFile=$P.bed --bimFile=$P.bim --famFile=$P.fam --AlleleOrder=alt-first"
  case "$wl" in
    binary)        model=binary_sparseGRM_a;  args="" ;;
    binary_sparse) model=binary_sparseGRM_vr; args="$SG --is_fastTest=FALSE" ;;
    quant)         model=quant_sparseGRM;     args="" ;;
    # The same markers through the PGEN reader (S8). Its workers each hold their
    # own pgenlib reader over shared metadata, which is a different sharing
    # arrangement from PLINK's one-mmap-many-readers, so it needs its own run.
    pgen)          model=binary_sparseGRM_a;  args=""
      bed="--pgenFile=$P.pgen --pvarFile=$P.pvar --psamFile=$P.psam" ;;
    pgen_quant)    model=quant_sparseGRM;     args=""
      bed="--pgenFile=$P.pgen --pvarFile=$P.pvar --psamFile=$P.psam" ;;
    # Controlled ER-eligible fraction (harness/make_rare_fixture.py). rare100 is
    # the workload that matters once the ER deferral is gone: EVERY marker takes
    # the exact test on a worker, where the bundled binary fixture puts only 4.1%
    # of markers there. A clean TSAN run on `binary` alone would be a statement
    # about 4% of the change.
    # The region (group-test) marker loop, threaded 2026-09-06: golden P2's
    # command. Its model lives in extdata/output, not $GOLD/step1, and the
    # marker-path flags below (--minMAC, --markers_per_chunk) are harmless to
    # it; bytes are compared against this workload's own 1-thread run.
    region)        model=example_binary_sparseGRM
      args="$SG --groupFile=input/group_new_snpid.txt --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous --maxMAF_in_groupTest=0.0001,0.001,0.01 --is_fastTest=TRUE --chrom=1"
      bed="--bedFile=./input/genotype_100markers.bed --bimFile=./input/genotype_100markers.bim --famFile=./input/genotype_100markers.fam --AlleleOrder=alt-first" ;;
    # LEDGER #172 (2026-09-29): the exact-by-weight routing runs inside the
    # per-cell prepares, on workers. The `region` workload above is
    # beta-weighted and never enters the routed body, so a clean run on it says
    # nothing about that code. These two do: golden P69's command (FlexRV
    # BURDEN, 24 routed cells through the BURDEN prepare) and golden P71's
    # (custom weights through SKAT-O, the ultra-rare prepare).
    region_atom)   model=example_binary_sparseGRM
      args="$SG --groupFile=input/group_new_snpid_flexrv_atom.txt --flexRV_score=REVEL --flexRV_maxMAF=0.01 --annotation_in_groupTest=lof:missense --r.corr=1 --is_fastTest=FALSE --chrom=1"
      bed="--bedFile=./input/genotype_100markers.bed --bimFile=./input/genotype_100markers.bim --famFile=./input/genotype_100markers.fam --AlleleOrder=alt-first" ;;
    region_atom_skato) model=example_binary_sparseGRM
      args="$SG --groupFile=input/group_new_snpid_withWeights_atom.txt --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous --maxMAF_in_groupTest=0.0001,0.001,0.01 --is_fastTest=FALSE --chrom=1"
      bed="--bedFile=./input/genotype_100markers.bed --bimFile=./input/genotype_100markers.bim --famFile=./input/genotype_100markers.fam --AlleleOrder=alt-first" ;;
    rare[0-9][0-9][0-9])
      model=binary_sparseGRM_a; args=""
      local rp="$REPO/harness/work/rarefix/$wl"
      [ -f "$rp.bed" ] || { echo "missing $rp.bed -- run harness/make_rare_fixture.py" >&2; exit 1; }
      bed="--bedFile=$rp.bed --bimFile=$rp.bim --famFile=$rp.fam --AlleleOrder=alt-first" ;;
  esac
  local base="$RUNS/$wl.t$nt" include="" modeldir="$GOLD/step1"
  case "$wl" in region*) modeldir="$REPO/extdata/output" ;; esac
  if [ -n "$MAXMK" ]; then
    # First N variant IDs of this workload's fixture. pvar has a #CHROM header
    # and ID in column 3; bim has no header and ID in column 2.
    local inc="$RUNS/$wl.include.txt" src
    if [ ! -f "$inc" ]; then
      # awk-with-a-counter rather than `| head -n`: head closing the pipe sends
      # SIGPIPE to awk, and `set -o pipefail` turns that into a fatal 141.
      src=$(echo "$bed" | tr ' ' '\n' | sed -n 's/^--pgenFile=//p')
      if [ -n "$src" ]; then
        awk -v n="$MAXMK" '!/^#/ {print $3; if (++i >= n) exit}' "${src%.pgen}.pvar" > "$inc"
      else
        src=$(echo "$bed" | tr ' ' '\n' | sed -n 's/^--bedFile=//p')
        awk -v n="$MAXMK" '{print $2; if (++i >= n) exit}' "${src%.bed}.bim" > "$inc"
      fi
    fi
    include="--idstoIncludeFile=$inc"
  fi
  # `env` rather than a bare VAR=... prefix: these are arguments to a shell
  # FUNCTION, and a prefix on a function call leaks into the shell afterwards.
  run_with_timeout "$LIMIT" "$base.log" env \
    TSAN_OPTIONS="halt_on_error=0 history_size=7 exitcode=0 log_path=$base.tsan" \
    RCPP_PARALLEL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    pixi run --manifest-path "$REPO/pixi.toml" bash "$WORK/run_r.sh" -f step2_SPAtests.R --args \
    $bed --minMAF=0 --minMAC=1 --LOCO=FALSE --markers_per_chunk=10000 --nThreads=$nt $include $args \
    --GMMATmodelFile="$modeldir/$model.rda" \
    --varianceRatioFile="$modeldir/$model.varianceRatio.txt" \
    --SAIGEOutputFile="$base.txt"
  local st="$RWT_STATUS" races threaded
  races=$(cat "$base".tsan.[0-9]* 2>/dev/null | grep -c "WARNING: ThreadSanitizer: data race" || true)
  threaded=$(grep -c "in-process threads" "$base.log" || true)
  # A dead or hung arm is its own outcome. Its output file is partial, so the
  # byte comparison would be meaningless; say what happened instead of
  # comparing. Under `plant` this is the defect being as bad as it is licensed
  # to be; under `clean` it is a failure however many races were reported.
  if [ "$st" -ne 0 ]; then
    local how="died (exit $st)"
    [ "$st" -eq 124 ] && how="TIMED OUT after ${LIMIT}s"
    printf '  %-14s t=%-2s races=%-4s %s\n' "$wl" "$nt" "$races" "$how"
    # Only under `clean` is this a failure of the script: a plant is licensed to
    # crash, and that is the whole point of running one.
    # An explicit `if` and `return 0`, not `[ clean ] && DIED=1; return`: under
    # `set -e` that `&&` list has status 1 whenever MODE is a plant, a bare
    # `return` propagates it, and errexit killed the script here after the
    # FIRST arm that died -- in the one mode where dying is the expected
    # outcome. The operator saw one arm line and exit 1 (2026-09-14 audit).
    if [ "$MODE" = clean ]; then DIED=1; fi
    return 0
  fi
  # A threaded arm that quietly ran serially reports zero races and proves
  # nothing -- the LEDGER #30 shape, and the reason the standing test asserts
  # this same line.
  if [ "$nt" -gt 1 ] && [ "$threaded" -eq 0 ]; then
    echo "  $wl t=$nt: DID NOT THREAD (no thread-count line) -- result is meaningless"
    return 0
  fi
  local ident="n/a"
  if [ "$nt" -gt 1 ]; then
    if cmp -s "$RUNS/$wl.t1.txt" "$base.txt"; then ident="identical"; else ident="DIFFERS"; fi
  fi
  printf '  %-14s t=%-2s races=%-4s bytes-vs-serial=%s\n' "$wl" "$nt" "$races" "$ident"
}

DIED=0
echo "=== TSAN runs ($MODE); ${MAXMK:-128,868} markers each"
for wl in ${SAIGE_TSAN_WORKLOADS:-binary binary_sparse quant region}; do
  for nt in 1 4 8; do tsan_run "$wl" "$nt"; done
done
echo
if [ "$MODE" = clean ]; then
  echo "Expected: races=0 everywhere, bytes-vs-serial=identical, and no arm died."
elif [ "$MODE" = plantreader ]; then
  echo "Expected (on SAIGE_TSAN_WORKLOADS=pgen): races>0 at t>1 and bytes-vs-serial=DIFFERS."
  echo "If this mode reports races=0, the CLEAN pgen result above gates nothing."
else
  echo "Expected: races>0 at t>1 (worker-vs-worker) and bytes-vs-serial=DIFFERS."
  echo "If this mode reports races=0, the CLEAN result above gates nothing."
  echo "An arm that died or TIMED OUT is also the plant working -- a shared"
  echo "MarkerScratch may crash before it finishes reporting."
fi
echo "reports and outputs: $RUNS"
exit $DIED
