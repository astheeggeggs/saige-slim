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

## 5. Investigating a saddlepoint failure on the real data (image a6213623 or later)

The sidecar names the step that failed; the dump names the vector. With
`SAIGE_SPA_DUMP_DIR` set, every failed saddlepoint writes
`<dir>/spaFail_<k>.tsv` (one row per sample: gtilde, mu, y; header: q, qinv,
m1, gneg, gpos, Tstat, var1, reason) and the sidecar's Note ends in
`dump=<file>`. Individual-level data: it stays in the workspace.

Wrapper: two lines. At the top, `SPA_DUMP="${SPA_DUMP:-}"` and
`OUT_FILE_SPADUMP="${OUT_FILE_SPADUMP:-}"`; before the run,
`if [[ -n ${SPA_DUMP} ]]; then mkdir -p spadump; export SAIGE_SPA_DUMP_DIR=$PWD/spadump; fi`;
after the deliveries,
`if [[ -n ${OUT_FILE_SPADUMP} ]]; then tar czf spadump.tgz spadump 2>/dev/null || tar czf spadump.tgz --files-from /dev/null; deliver spadump.tgz "${OUT_FILE_SPADUMP}"; fi`.
`submit()`: `"--output", f"OUT_FILE_SPADUMP={od}/{pre}.spaDump.tgz"` and
`env={"SPA_DUMP": "1", ...}`. Run ONE trait (the 24-fallback one) with
`REGIONS` set to the failing regions only.

Analysis runs in a dsub job on the SAIGE image (the notebook VM has no SAIGE;
the exported kernels are what reproduce the finder exactly):

```r
# Rscript inside the SAIGE image, one dump file as input; prints summary numbers only
suppressMessages(library(SAIGE))
f <- commandArgs(TRUE)[1]
hdr <- strsplit(readLines(f, n = 1), "\t")[[1]]; d <- read.table(f, header = TRUE, skip = 1, sep = "\t")
q <- as.numeric(hdr[2]); qinv <- as.numeric(hdr[4]); tol <- .Machine$double.eps^0.25
mu <- d$mu; g <- d$gtilde
r1 <- getroot_K1_Binom(0, mu, g, q, tol, 1000L); r2 <- getroot_K1_Binom(0, mu, g, qinv, tol, 1000L)
cat("N", length(g), "carriers", sum(g > 0.5 * max(g)), "gneg", sum(g[g < 0]), "gpos", sum(g[g > 0]), "qinv", qinv, "\n")
cat("root q:", r1$root, r1$Isconverge, r1$niter, " root qinv:", r2$root, r2$Isconverge, r2$niter, "\n")
for (t in -2^(0:24)) cat(sprintf("t=%-10.3g K1-qinv=%-12.4g K2=%-10.3g\n", t, K1_adj_Binom(t, mu, g, qinv), K2_Binom(t, mu, g)))
p <- SPA_pval(mu, g, q, qinv, 1e-100, tol, FALSE, "binary", FALSE); cat("composed p:", p, "\n")
```

What to look for: where K1 - qinv changes sign (the root's location against
the bracket cap 2^20 / max|g|), whether K2 underflows to 0 before it, and
whether any entry makes K1 non-finite. That names the property my synthetic
vectors lack.

**What it found (2026-09-27, `aou_spa_dumpcheck.sh` on the 20 BenIntNeo dumps,
LEDGER #169):** `root(q)` converged in every dump, the mirror was unattainable
(0 iterations), and `Korg(root1)` was Inf with `K2(root1)` finite -- the CGF's
`exp(g t)` overflows past ~709 while K1/K2 are written in `exp(-g t)`. Fixed
in the kernels (commit `4ff9576b`, image `1.5.2-dev-4ff9576b` and later). The line to
read in this log is the `q-side saddle:` one: `isSaddle=FALSE` with
`Korg(root1)=Inf` is #169; `root(qinv)` with `niter=0` is the #168 (a) rule,
not the finder. The sidecar's Note on that image and later carries
`qsaddle=INVALID` directly.

## 6. Rolling the fix across the ancestries (2026-09-27, image `1.5.2-dev-3e92d89d`)

**Use `1.5.2-dev-3e92d89d` or later, not `4ff9576b`**: on `4ff9576b` a chromosome whose rerun
list holds only regions without a pLoF variant dies on the `--flexRV_lofAnno`
file-level check (LEDGER #170) while the job reports SUCCESS -- sas/BenIntNeo
lost chromosome 5 that way on both `8ca5d0e4` and `4ff9576b`. The stage
checker reads `chrStatus` and flags it.

Order decided by the owner: **sas/BenIntNeo first** (the one trait still
unspliced, and the proof the #169 image lands: 24 fallbacks on `8ca5d0e4`,
expected 0), then **eas and afr** together, then **amr**, then **eur**. Only
the FlexRV arm of the BINARY traits reruns (a quantitative trait has no
saddlepoint; the SKAT-O arms scanned clean everywhere); each job is filtered
to its scan's regions (`regions_for`), so eur's size is the number of failing
regions, not the cohort. amr: copy v4's flags exactly, which for amr means no
`--sparseGRMFile` (its GRM is not sparse; the step-2 solve does not finish).
eur: any job expected past ~45 min goes on demand, not Spot (preemption is
correlated by zone).

```python
IMG7 = "gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-3e92d89d"
for anc in ["sas"]:                      # then ["eas", "afr"], then ["amr"], then ["eur"]
    for pre, pheno in binary_traits(anc):        # the v4 job prefixes with a non-empty regions list
        submit_rerun(anc, pheno, pre, STAMP[anc], IMG7)
```

**Per-stage check, before the next ancestry is launched.** "Every
saddlepoint converged, or is NA for a stated reason" is read off the
delivered `<pre>.spaFallbacks.txt.gz` sidecars (the wrapper gzips them; a 0-fallback run delivers a header-only placeholder of ~92 bytes, and the checker reads `.txt` or `.txt.gz`), which are the COMPLETE record: a
region-path cell whose p is NA is not written to the results file at all
(LEDGER #164; the splice's `ROWS` line makes the missing rows visible), so
the results file cannot show it and the sidecar is where every failure is
named. Copy the stage's sidecars to a local directory and run

```sh
bash harness/tools/aou_spa_stagecheck.sh <dir>      # one line per trait, one verdict
```

which prints per trait `fallbacks F (exact p E, NA X)` with the histogram of
`spa: <reason>` from the Note column, lists every NA row's Note, and ends in
`STAGE OK` (no NA) or `STAGE: N NA cell(s) -- read their Notes` (exit 1). On
this image the expectation is zero fallbacks per trait, and a fallback that
ran the exact test is a valid p (it is counted, not flagged). An NA is
acceptable only when its Note says why (`carriers, above the fallback's
exact-test ceiling` after a reason other than #169's `saddle(q) invalid`
with `qsaddle=INVALID`); a `qsaddle=INVALID` on this image is a new
mechanism and stops the rollout. A `CHROMOSOME FAILED` line is a lost
chromosome: its regions are absent from the results and the trait is not
spliceable until rerun; the job log's SAIGE output for that chromosome says why. Then splice (section 3; every `changed
row(s), all explained` line green, no `WARNING ... came back identical`),
and only then launch the next ancestry.

**Positive control (afr/BenIntNeo).** The genes the scan reported at
p ~1e-174 were the #163/#164 normal-p artefact; on this image they should
come back at ~1e-3 (the dominant carrier's fitted probability times the
rest's tail). If they come back at 1e-100 or NA, stop.
