# Filtered rerun after LEDGER #163/#164 (AoU, 2026-09-26)

The scan (`aou_spa_scan.sh`) leaves, per delivered FlexRV file,
`rerun/<pre>.regions.txt`: the regions with a failed saddlepoint plus the
regions with any burden cell below the tail threshold. Only those regions are
rerun, on the fixed image, and their rows are spliced over the old file.
Sound because the fix changes ONLY tests whose saddlepoint failed: every other
row comes back byte-identical, which the splice checks.

## 1. Wrapper patch (saige_step2_brava_main_original.sh)

A `REGIONS` env var (comma-separated region IDs) restricts every
chromosome's group file to those regions before SAIGE runs; a chromosome
with none of them is skipped, and the skip is recorded in chrStatus.
Apply after the block that resolves `CHR_GRP[$c]` (the loop
`for c in "${CHR_LIST[@]}"; do ... done` that validates every chromosome):

```bash
# ---- LEDGER #164 filtered rerun: keep only the listed regions --------------
REGIONS="${REGIONS:-}"                  # "ENSG...,ENSG...": rerun these regions only
SKIPPED_CHRS=()
if [[ -n ${REGIONS} && ${TESTTYPE} == "group" ]]; then
  tr ',' '\n' <<< "${REGIONS}" | sed '/^$/d' | sort -u > regions.keep
  echo "REGIONS: $(wc -l < regions.keep) region(s) to rerun"
  kept=()
  for c in "${CHR_LIST[@]}"; do
    sub="group.chr${c}.subset.txt"
    gzip -cdf "${CHR_GRP[$c]}" | awk 'NR==FNR{r[$1]=1; next} ($1 in r)' regions.keep - > "${sub}"
    if [[ -s ${sub} ]]; then
      CHR_GRP[$c]="${sub}"; kept+=("${c}")
      echo "chr${c}: $(awk '$2=="var"' "${sub}" | wc -l) region(s) kept"
    else
      SKIPPED_CHRS+=("${c}")
    fi
  done
  (( ${#kept[@]} > 0 )) || die "REGIONS matched no region in any chromosome's group file"
  CHR_LIST=("${kept[@]}")
  (( ${#CHR_LIST[@]} > 1 )) && MULTI_CHR=1 || MULTI_CHR=""
fi
```

and add `skipped_no_region\t${SKIPPED_CHRS[*]:-none}` to the chrStatus
printf near the end. The FlexRV score-line check that follows runs on the
subset, so it still holds one score line per region.

## 2. Launch (notebook), one job per trait with a non-empty regions list

```python
def regions_for(anc, arm_dir, pre, stamp):
    p = f"{BUCKET}/SAIGE_SLIM_TESTING/spa_scan/{anc}/{arm_dir}/{stamp}/rerun/{pre}.regions.txt"
    txt = subprocess.run(["gcloud", "storage", "cat", p], capture_output=True, text=True).stdout
    return [r for r in txt.split() if r]

def submit_rerun(anc, pheno, pre, stamp, image, trait="binary", **kw):
    regs = regions_for(anc, "flexrv_v4", pre, stamp)
    if not regs:
        print(pre, "nothing to rerun"); return None
    print(pre, len(regs), "regions")
    return submit(name=f"fxv5-{anc}-{tag(pheno)}", image=image,
                  masks="pLoF:damaging_missense:other_missense", fasttest="auto",
                  outdir=f"{anc}/flexrv_v5", anc=anc, pheno=pheno, trait=trait,
                  chrs=range(1, 23), plk=lambda c, a=anc: cod(c, a), grp=am_grp,
                  pooled=True, sidecars=True,
                  env={"FLEXRV_SCORE": "AM", "FLEXRV_MAXMAF": "0.001", "MAXMAF": "0.001",
                       "RCORR": "1", "REGIONS": ",".join(regs)}, **kw)
```

Use exactly the v4 job's masks, MAXMAF, model and cutoff (copy them from
`RUNS["fxv4-..."]["masks"]`, `["env"]`) -- the only differences must be the
image, `REGIONS`, and the output directory. `fasttest="auto"` is what v4
used (group runs get no --is_fastTest).

## 3. Splice + controls (in-cloud, `aou_spa_splice.sh`)

Replaces the rerun regions' rows in the old results file, singleAssoc and
pooledTests with the new ones, and checks two things before writing:
- CONTROL: a rerun region with no failing row in `failures.tsv` must come
  back byte-identical (every row), or the splice refuses -- that is the
  proof that only failed tests moved;
- every region with a failing row must have at least one changed row.
Writes `<pre>.txt.gz` etc. under `<anc>/flexrv_v4_fixed/` plus
`<pre>.splice.txt` naming the regions replaced and the row counts.

## 4. Per-chromosome timeout (diagnostic for a chromosome that never returns)

On the first SAS rerun, HTN's chr16 (5 regions) took 12 minutes where every
other chromosome took 5-20 s (27 fallbacks, one cell with 209 carriers in the
exact test's resampling regime; it finished). The data property behind the
slowdown is not identified: the
exact test's cost is not set by carriers or weights (measured at N = 50,000,
13-209 carriers, weights 1 or 2e-6, all under 5 s). The wrapper prints a
chromosome's SAIGE log only on failure, so a hang leaves nothing. A timeout
turns it into a failure whose log ends with the LAST `Analyzing Region` line,
which names the region. In `launch_chr`, wrap the command:

```bash
  launch_chr() {                 # $1 = chromosome; sets LAST_PID
    local c="$1" cmd_c
    cmd_c=("${cmd[@]}" ${CHR_GENO[$c]} --chrom="chr${c}")
    [[ ${TESTTYPE} == "group" ]] && cmd_c+=(--groupFile="${CHR_GRP[$c]}")
    if [[ -n ${CHR_TIMEOUT_SECONDS:-} ]]; then
      timeout --signal=TERM --kill-after=30 "${CHR_TIMEOUT_SECONDS}" \
        "${cmd_c[@]}" --SAIGEOutputFile="${OUTBASE}.chr${c}" > "saige.chr${c}.log" 2>&1 &
    else
      "${cmd_c[@]}" --SAIGEOutputFile="${OUTBASE}.chr${c}" > "saige.chr${c}.log" 2>&1 &
    fi
    LAST_PID=$!
  }
```

(`timeout` is GNU coreutils and is in the image.) Then rerun HTN's chr16
alone with a short limit and the same region list:

```python
r = RUNS["fxv4-sas-htn"]; env = dict(r["env"]); env["CHR_TIMEOUT_SECONDS"] = "600"
env["REGIONS"] = ",".join(regions_for("sas", "fxv4-sas-htn.HTN.sas.chr1-22", STAMP))
submit(name="diag-sas-htn-chr16", image=IMG5, masks=r["masks"], fasttest=r["fasttest"],
       outdir="sas/diag", anc="sas", pheno="HTN", trait="binary", firth=r["firth"], env=env,
       geno=r["geno"], chrn="16", plk=cod(16, "sas"), grp=am_grp(16), testtype="group",
       pooled=False, sidecars=True, sparse_grm=r["sparse_grm"])
```

The job fails at 10 minutes and its log then holds `saige.chr16.log` in full;
`grep "Analyzing Region" | tail -1` names the region, and the same region's
144 v4 rows (MAC, Number_ultra_rare, Number_rare per set) are the next thing
to read. A region-level p is aggregate output.
