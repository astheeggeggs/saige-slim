#!/bin/bash
# How many samples of each null model are PINNED (fitted probability at the 0/1
# boundary), read from the step-2 logs already delivered -- no model is opened,
# no new SAIGE run (LEDGER #174). Since the fit gates (R/FitGates.R) every
# step-2 run prints, when any fitted value is below 1e-8 or above 1 - 1e-8,
#   the null model in <file> has K of N fitted values (x%) at the 0/1 boundary (range a .. b)
# once per chromosome. This job greps the newest log of every job whose name
# starts with PREFIX and writes ONE table, <OUT_DIR>/pinned.<tag>.tsv:
#   anc  trait  job  status  pinned  N  percent  min_fitted  max_fitted
# status = warning (the line was found), none_reported (no such line: no pinned
# sample at 1e-8, or a run with the gates off / an image older than the gates --
# `gates` says which the log suggests), no_log.
# Runs INSIDE a dsub job (gcr.io/google.com/cloudsdktool/google-cloud-cli:slim).
# A count from 1 to MASK-1 is written as "<MASK" and its percentage as "masked":
# the table is safe to read and to keep. Env:
#   LOGS     gs://<bucket>/dsub/logs   (or a local dir of the same shape, for the self-test)
#   PREFIX   job-name prefix, default fxv4-   (names are <prefix><anc>-<trait>[-b<bin>])
#   TAG      label for the output file, default v4
#   OUT_DIR  local dir dsub delocalizes (--output-recursive)
#   MASK     default 20
set -euo pipefail
: "${LOGS:?}"; : "${OUT_DIR:?}"; PREFIX="${PREFIX:-fxv4-}"; TAG="${TAG:-v4}"; MASK="${MASK:-20}"
mkdir -p "${OUT_DIR}"
if [[ ${LOGS} == gs://* ]]; then
  lsall(){ gcloud storage ls "${LOGS}/${PREFIX}*/**/*.log" 2>/dev/null || true; }
  catf(){ gcloud storage cat "$1"; }
else
  lsall(){ find "${LOGS}" -type f -name '*.log' -path "${LOGS}/${PREFIX}*" 2>/dev/null || true; }
  catf(){ cat "$1"; }
fi
lsall | sort -r > logs.list          # newest date/time directory of a job first
echo "logs listed: $(( $(wc -l < logs.list) ))"
# one row per job name: the path component right after LOGS/
awk -v base="${LOGS}/" '{p=$0; sub(base, "", p); n=split(p, a, "/"); print a[1] "\t" $0}' logs.list > joblogs.tsv
cut -f1 joblogs.tsv | sort -u > jobs.list
echo "jobs: $(( $(wc -l < jobs.list) ))"
out="${OUT_DIR}/pinned.${TAG}.tsv"
printf 'anc\ttrait\tjob\tstatus\tgates\tpinned\tN\tpercent\tmin_fitted\tmax_fitted\n' > "${out}"
while read -r job; do
  [[ -n ${job} ]] || continue
  rest="${job#"${PREFIX}"}"; anc="${rest%%-*}"; trait="${rest#*-}"; trait="${trait%-b[0-9]*}"
  line=""; gates="on_or_unknown"; seen=0
  while IFS=$'\t' read -r j f; do
    [[ ${j} == "${job}" ]] || continue
    seen=1
    # the pattern, then (only if absent) whether the log says the gates were off
    line=$(catf "${f}" 2>/dev/null | grep -m1 -E 'fitted values \([^)]*\) at the 0/1 boundary' || true)
    if [[ -n ${line} ]]; then break; fi
    if catf "${f}" 2>/dev/null | grep -m1 -q -E 'SAIGE_FIT_GATES=0|fit gates are OFF|gates are off'; then gates="off"; fi
  done < joblogs.tsv
  if (( ! seen )); then
    printf '%s\t%s\t%s\tno_log\tNA\tNA\tNA\tNA\tNA\tNA\n' "${anc}" "${trait}" "${job}" >> "${out}"; continue
  fi
  if [[ -z ${line} ]]; then
    printf '%s\t%s\t%s\tnone_reported\t%s\t0\tNA\tNA\tNA\tNA\n' "${anc}" "${trait}" "${job}" "${gates}" >> "${out}"; continue
  fi
  echo "${line}" | sed -E 's/.*has ([0-9]+) of ([0-9]+) fitted values \(([^%]+)%\) at the 0\/1 boundary \(range ([^ ]+) \.\. ([^) ]+)\).*/\1\t\2\t\3\t\4\t\5/' | \
    awk -F'\t' -v OFS='\t' -v anc="${anc}" -v tr="${trait}" -v job="${job}" -v mask="${MASK}" '
      { k=$1+0; shown=k; pct=$3; if(k > 0 && k < mask){ shown="<" mask; pct="masked" }
        print anc, tr, job, "warning", "on", shown, $2, pct, $4, $5 }' >> "${out}"
done < jobs.list
# one row per (anc, trait): the bins of a trait share its model
{ head -1 "${out}"; tail -n +2 "${out}" | sort -t$'\t' -k1,1 -k2,2 -k4,4r | awk -F'\t' '!s[$1 FS $2]++'; } > "${out}.dedup" && mv "${out}.dedup" "${out}"
n=$(($(wc -l < "${out}") - 1))
# counts through awk's END, not `wc -l`: BSD wc pads its output
nw=$(awk -F'\t' 'NR>1 && $4=="warning"{c++} END{print c+0}' "${out}"); nm=$(awk -F'\t' -v m="<${MASK}" 'NR>1 && $6==m{c++} END{print c+0}' "${out}")
nn=$(awk -F'\t' 'NR>1 && $4=="none_reported"{c++} END{print c+0}' "${out}")
echo "models: ${n}; with pinned samples: ${nw} (of which masked below ${MASK}: ${nm}); none reported: ${nn}"
gzip -f "${out}"
echo done
