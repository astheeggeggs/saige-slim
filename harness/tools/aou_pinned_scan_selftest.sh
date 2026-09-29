#!/bin/bash
# Self-test for aou_pinned_scan.sh on FABRICATED logs (no AoU data): the warning
# is found and parsed, a count below the mask is written "<20" with its
# percentage masked, a log without the line is none_reported, the newest log of
# a resubmitted job wins, the bins of one trait collapse to one row, and a
# gates-off log is flagged.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d); trap 'rm -rf "${T}"' EXIT
cd "${T}"
mk(){ mkdir -p "$(dirname "$1")"; printf '%s\n' "${@:2}" > "$1"; }
W(){ echo "Warning message: the null model in /x/$1.rda has $2 of $3 fitted values ($4%) at the 0/1 boundary (range $5 .. $6). That is what separation in the covariate model does"; }
mk logs/fxv4-xx-traita/u/20260901/100000/j-task-None.log "start" "$(W traita 77 4089 1.88 2.22e-16 0.41)" "end"
mk logs/fxv4-xx-traitb/u/20260901/100000/j-task-None.log "start" "$(W traitb 5 9000 0.0556 3e-12 0.3)" "end"
mk logs/fxv4-xx-traitc/u/20260901/100000/j-task-None.log "start" "Analysis done" "end"
mk logs/fxv4-xx-traitd/u/20260901/100000/j-task-None.log "start" "$(W traitd 900 9000 10 1e-16 0.9)" "end"      # an OLD submission
mk logs/fxv4-xx-traitd/u/20260902/110000/j-task-None.log "start" "$(W traitd 30 9000 0.333 1e-16 0.9)" "end"    # the newest wins
mk logs/fxv4-yy-traite-b0/u/20260901/100000/j-task-None.log "start" "$(W traite 40 50000 0.08 1e-15 0.5)" "end"
mk logs/fxv4-yy-traite-b1/u/20260901/100000/j-task-None.log "start" "$(W traite 40 50000 0.08 1e-15 0.5)" "end"
mk logs/fxv4-yy-traitf/u/20260901/100000/j-task-None.log "start" "WARNING: SAIGE_FIT_GATES=0 -- the fit gates are OFF" "end"
mk logs/skv4-xx-traita/u/20260901/100000/j-task-None.log "start" "$(W traita 77 4089 1.88 2.22e-16 0.41)"       # another prefix: not read
LOGS="${T}/logs" OUT_DIR="${T}/out" bash "${HERE}/aou_pinned_scan.sh" > run.log 2>&1 || { cat run.log; echo "SELFTEST FAILED: script exited non-zero"; exit 1; }
gzip -dc out/pinned.v4.tsv.gz > out.tsv
fail=0
check(){ if awk -F'\t' "$1" out.tsv | grep -q .; then echo "ok   $2"; else echo "FAIL $2"; fail=1; fi; }
check 'NR>1 && $1=="xx" && $2=="traita" && $4=="warning" && $6=="77" && $7=="4089" && $8=="1.88" && $9=="2.22e-16"' "a count at or above the mask is written with its percentage and range"
check 'NR>1 && $2=="traitb" && $6=="<20" && $8=="masked" && $7=="9000"' "a count below the mask is masked, and so is its percentage"
check 'NR>1 && $2=="traitc" && $4=="none_reported" && $5=="on_or_unknown" && $6=="0"' "no warning line: none_reported"
check 'NR>1 && $2=="traitd" && $6=="30"' "the newest log of a resubmitted job wins"
check 'NR>1 && $1=="yy" && $2=="traite" && $6=="40"' "a binned trait is parsed"
[[ $(awk -F'\t' 'NR>1 && $2=="traite"' out.tsv | wc -l) -eq 1 ]] && echo "ok   the bins of one trait are one row" || { echo "FAIL bins not collapsed"; fail=1; }
check 'NR>1 && $2=="traitf" && $4=="none_reported" && $5=="off"' "a gates-off log is flagged"
[[ $(($(wc -l < out.tsv) - 1)) -eq 6 ]] && echo "ok   6 models (the other prefix is not read)" || { echo "FAIL row count $(wc -l < out.tsv)"; fail=1; }
grep -q "models: 6; with pinned samples: 4 (of which masked below 20: 1); none reported: 2" run.log && echo "ok   summary line" || { echo "FAIL summary line: $(tail -3 run.log | head -1)"; fail=1; }
! grep -q $'\t5\t9000' out.tsv && echo "ok   the masked count appears nowhere in the table" || { echo "FAIL masked count leaked"; fail=1; }
[[ ${fail} -eq 0 ]] && echo "selftest OK: 10 checks" || { echo "SELFTEST FAILED"; cat out.tsv; exit 1; }
