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

Use exactly the v4 job's masks, MAXMAF, model, cutoff **and `geno`** (copy them
from `RUNS["fxv4-..."]`) -- the only differences must be the image, `REGIONS`,
and the output directory. **The v4 records carry `fasttest` `FALSE`, not `auto`**
(read from `RUNS["fxv4-<anc>-<trait>-b0"]` for amr and eur, 2026-10-02; the launch
set `FASTTEST="FALSE"`, `docs/state/aou-saige-parameters.md`). On a GROUP run the two
are the same command: the wrapper's `auto` gives `--is_fastTest=TRUE` to variant
tests only, so `fasttest="auto"` in the cells here passes no flag. Prefer
`fasttest=r["fasttest"]` from the record, as the diagnostic cell in section 4 does:
it stays right if a cell is ever reused for a variant run, where `auto` and `FALSE`
differ.

**`submit()` DEFAULTS `geno` TO `bed`, AND EVERY PRODUCTION JOB USES PGEN.**
Measured 2026-10-01, the hard way: a cell that did not pass `geno` ran for 61 s
and exited 1, because the coding subsets exist only as
`...cod.{pgen,pvar,psam}` and there is no `.bed/.bim/.fam` at that prefix. The
failure is NOT distinctive -- Batch reports "exit code 1" and the duration looks
like a localization problem -- so the cheapest check on any new submit cell is
`RUNS[name]["geno"]` after submitting, which should read `pgen`.

PGEN is also what you want on cost grounds, and for a reason that is not decode
speed: a dsub worker LOCALIZES THE WHOLE GENOTYPE FILE before the script starts,
and the AoU coding subsets are ~180x smaller as PGEN (2.3 GiB -> 13 MB on chr20
afr). The compute difference on a gene test is only 1.13x; the data movement is
the win. **Pass `geno="pgen"` explicitly in every new cell**, or better, change
the notebook's default -- the current one is a footgun that costs a job and ~60 s
to discover.

## 3. Splice + controls (in-cloud, `aou_spa_splice.sh`)

Replaces the rerun regions' rows in the old results file, singleAssoc and
pooledTests with the new ones, and checks two things before writing:
- CONTROL: a rerun region with no failing row in `failures.tsv` must come
  back byte-identical (every row), or the splice refuses -- that is the
  proof that only failed tests moved;
- every region with a failing row must have at least one changed row.
Writes `<pre>.txt.gz` etc. under `<anc>/flexrv_v4_fixed/` plus
`<pre>.splice.txt` naming the regions replaced and the row counts.

**The singleAssoc holds two kinds of row (found on the eas/afr stage,
2026-09-27).** A pseudo-marker row's MarkerID (`Region:Group:maxMAF:Weight`)
names its region; a REAL marker -- rare but not ultra-rare, MAC above the
collapse cutoff, which exists whenever `2 * N * maxMAF > 10`, i.e. N above
~5,000 at max MAF 0.001 (eas, afr, amr, eur; never sas) -- is tested singly
under its own `chr:pos:ref:alt` ID. The splice keys those by MarkerID: the
rerun's rows replace the old rows for the same markers, a changed one must be
a recorded real-marker failure or a marker at MAC <= `ER_MAC`, and the log
prints `real-marker rows: K rerun marker(s), ...`. Before this the splice
refused every eas/afr job with `ERROR: N new rows outside the region list`.
A real marker sits in every region that spans it and its row is written once
per region (the protocadherin cluster at chr5:140.98 Mb puts one variant in a
dozen), so the splice replaces as many old copies as the rerun wrote and keeps
the rest for the regions not rerun (`... Z old cop(y/ies) kept`); a rerun row
that differs from its copies is refused like any other unexplained change.
Also: a recorded failure is the pseudo-marker's OWN single test, so in the
results file a burden cell with `Number_rare > 0` is a different test that
need not have changed -- that is a `NOTE`, and the exact check (the
pseudo-marker row itself changed) is made in the singleAssoc. A `ROWS
.pooledTests` line is expected wherever the v4 artefact had fired the CCT
fallback rule: a plain Cauchy discloses nothing and writes no pooled row.

**A burden cell with rare, non-collapsed markers is a test the scan could not
observe** (eas/T2Diab, same evening): its vector is w*(real markers) +
pseudo-marker, its saddlepoint can fail while every component's own single
test converged, and v4 recorded nothing (pre-#164 the normal p was reported
silently). Seen as a MAC-47 `sbeta20` cell at 1.6e-7 in v4 and 0.068 on the
fixed image, same BETA. The splice accepts a changed cell with
`Number_rare > 0` at REGION level (the region is in the rerun list for a
failure or a tail p; the tail arm exists for exactly this), keeps the exact
cell-level rule for every observable cell, and prints `unobservable burden
cells (Number_rare > 0) accepted at region level: N (p rose in U, fell in D)`.
Read that line: the artefact p's were too SMALL, so U should dominate; a
stage where D is large is a question, not a pass.

**A pseudo-marker's single test can change without a recorded failure**
(eas/Asth, afr/BenIntNeo, same evening): since LEDGER #166 every set's
pseudo-marker takes the exact test when its OWN unweighted allele count is
<= `--max_MAC_for_ER`, where v4 gave that to `unw_unw` alone. That count is
neither the set's MAC in the results (which adds the rare, non-collapsed
markers) nor a weighted set's `AC_Allele2` (a weighted sum). The splice
reads it exactly from the set's 0/1-weight TWIN: a set with 0/1 weights
(`unw`|`lof` score x `unw`|`MACk` MAF) reports its own count in its
pseudo-marker's AC (LEDGER #133), and every set collapses the same markers
as -- or a subset of -- the 0/1 set with its zero pattern (`lof x mafbeta`
keeps the LoF markers, as `lof x unw` does; `cbrt x MAC25` keeps the MAC <= 25
markers, as `unw x MAC25` does), so the twin's count bounds it from above.
Twin at or below the cutoff explains the row, printed as `pseudo-marker rows
explained by their own exact-test eligibility (LEDGER #166, ...): K`; a
pseudo-marker changing where its twin's count is above the cutoff is still
refused. The `unw_unw` bound alone missed `lof x mafbeta` on three afr
traits, because the LoF subset can be at the cutoff while the whole
collapsed set is not.


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
regions, not the cohort. amr: copy v4's flags exactly -- it runs WITH a sparse
GRM like every other ancestry (the 2026-09-21 "drop `--sparseGRMFile` for amr"
prescription is void: the GRM was rebuilt with the relatedness cutoff changed,
owner, 2026-09-28; `submit_rerun` copies `sparse_grm` from the v4 record and
never needs to know).
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

## 7. The plotting pass: gene-level summaries first (2026-09-28)

Never plot from the per-trait results files (the 2026-09-26 egress warning
coincided with that). `harness/tools/aou_flexrv_summary.sh`, run as ONE dsub
job per ancestry on the cloud-cli image, streams every FlexRV and SKAT-O
results file of the ancestry and writes one small gzipped TSV to the bucket:
one row per (trait, gene), gene level only -- `p_flexrv` (the FlexRV file's
Cauchy row), `p_flexrv_best`/`set_best`/`n_sets`, `p_skato` (the SKAT-O
file's Cauchy row over its 18 cells), `p_skato_best`/`cell_best`,
`p_burden_best`, `p_skat_best`, and `source` = `fixed` when the trait's job
was spliced (`flexrv_v4_fixed/`) else `v4` (nothing to rerun, or a
quantitative trait). A gene in one arm only carries NA on the other. A second table per
ancestry, `<anc>.skato_cells.v4.tsv.gz`, holds every SKAT-O cell of every gene
wide (`burden:<mask>@<maxMAF>`, `skato:<mask>@<maxMAF>`), which the per-mask
QQ grids are drawn from (a few tens of MB gzipped for eur). Self-test:
`aou_flexrv_summary_selftest.sh`. The notebook then loads five files of a
few MB, and displays figures only.

```python
def submit_summary(anc, image="gcr.io/google.com/cloudsdktool/google-cloud-cli:slim"):
    args = ["aou_dsub", "--name", f"summary-{anc}", "--image", image,
            "--env", f"ANC={anc}", "--env", f"BASE={BUCKET}/SAIGE_SLIM_TESTING", "--env", "RUNTAG=v4",
            "--output-recursive", f"OUT_DIR={BUCKET}/SAIGE_SLIM_TESTING/summary_v4",
            "--machine-type", "n2-standard-2", "--disk-size", "20", "--script", "aou_flexrv_summary.sh"]
    script = "source /home/jupyter/aou_dsub.bash\n" + " ".join(shlex.quote(a) for a in args)
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=os.environ)
    print(anc, "rc =", out.returncode, "|", ((out.stdout or "") + (out.stderr or "")).strip()[-400:])
for anc in ANCS: submit_summary(anc)
```

**Rendering inside dsub instead of the notebook** (`aou_flexrv_qq.py`, image
`gcr.io/weighty-elf-452116-c7/pyplot:20260928` from `docker/Dockerfile.pyplot`):
the same grids, drawn on a batch VM from the two tables as dsub inputs, PNGs
delivered under `summary_v4/qq/<anc>/`; nothing is read into or displayed
from the notebook. Self-test `aou_flexrv_qq_selftest.sh <python>`.

```python
def submit_qq(anc, image="gcr.io/weighty-elf-452116-c7/pyplot:20260928"):
    sv = f"{BUCKET}/SAIGE_SLIM_TESTING/summary_v4"
    args = ["aou_dsub", "--name", f"qq-{anc}", "--image", image, "--env", f"ANC={anc}",
            "--input", f"IN_SUMMARY={sv}/{anc}.flexrv_vs_skato.v4.tsv.gz", "--input", f"IN_CELLS={sv}/{anc}.skato_cells.v4.tsv.gz",
            "--output-recursive", f"OUT_DIR={sv}/qq/{anc}",
            "--machine-type", "n2-standard-2", "--disk-size", "20", "--script", "aou_flexrv_qq.py"]
    script = "source /home/jupyter/aou_dsub.bash\n" + " ".join(shlex.quote(a) for a in args)
    r = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=os.environ)
    print(r.stdout[-300:] or r.stderr[-500:])
```

### 7a. The old-vs-new comparison, per ancestry (2026-10-01)

Measures what LEDGER #172's fix actually moved. Run on eas 2026-09-30 and the
most informative read of that day; sas, amr, eur and afr are still owed one.
The handoff called this `submit_summary_prefix`, which is not a function that
exists -- it is `submit_summary` run a SECOND time with two env vars changed.

**How the two sides are produced.** `aou_flexrv_summary.sh` tags each row's
`source`: `v7` when the trait's file came from `flexrv_v7_checked/`, `fixed`
when it came from `flexrv_v4_fixed/`, `v4` when nothing was rerun (every
quantitative trait). Pointing `FLEXRV_V7_DIR` at a directory that does NOT
exist empties `fx7.list` (`:37`, `:42`) and the script falls back to the fixed
files -- so the same script produces the OLD side with no special mode.

```python
def submit_summary_prefix(anc, image="gcr.io/google.com/cloudsdktool/google-cloud-cli:slim"):
    args = ["aou_dsub", "--name", f"sumpre-{anc}", "--image", image,
            "--env", f"ANC={anc}", "--env", f"BASE={BUCKET}/SAIGE_SLIM_TESTING", "--env", "RUNTAG=v4",
            "--env", "FLEXRV_V7_DIR=__does_not_exist__",      # forces the fallback to flexrv_v4_fixed/
            "--output-recursive", f"OUT_DIR={BUCKET}/SAIGE_SLIM_TESTING/summary_v4_prefix",
            "--machine-type", "n2-standard-2", "--disk-size", "20", "--script", "aou_flexrv_summary.sh"]
    script = "source /home/jupyter/aou_dsub.bash\n" + " ".join(shlex.quote(a) for a in args)
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=os.environ)
    print(anc, "rc =", out.returncode, "|", ((out.stdout or "") + (out.stderr or "")).strip()[-400:])
for anc in ["sas", "amr", "eur", "afr"]: submit_summary_prefix(anc)
```

**A submit cell MUST print stderr and the return code.** `submit_summary`
printed `.stdout` alone for its whole life, so a submit that never launched
anything printed an empty-looking result and read as success -- observed
2026-10-01, when all four prefix jobs printed `[]` and the cause was on stderr.
`submit_check_v7` had it right (`out.stdout or out.stderr`); the two summary
cells did not, and now do.

**Read the job log's `flexrv: N v4 files, M fixed, K v7` line before using the
output.** It is printed by `aou_flexrv_summary.sh:43` to STDOUT, so it is in the
job's log file, not in the output TSV. The provider on this project is
**`google-batch`** (`google-cls-v2` is rejected by the installed dsub, which
offers only `local`, `google-batch`, `test-fails`):

```sh
dstat --provider google-batch --project "$GOOGLE_PROJECT" --location us-central1 \
      --names sumpre-<anc> --full          # the record, including its `logging` path
gcloud storage cat <logging-path>/sumpre-<anc>.log | grep '^flexrv:'
```

The merge's `assert (j.source_old == "fixed").any()` catches the same failure
from the other end, at the point where it would mislead, so the log line is the
cheaper check rather than the only one. On a PREFIX run `K` must be **0** -- a non-zero K means the directory
name resolved after all and you have built the new side twice. On the NEW side
`K` must equal the ancestry's job count (14c). This is the only thing that
distinguishes the two tables: same filename, same `RUNTAG`, different directory.

**The comparison, with the assert that makes it mean anything.** A
self-comparison gives a median ratio of exactly 1.0000 -- and so does the
correct comparison, because the fix is surgical. The number cannot tell them
apart, so assert on the `source` columns first and never skip it.

```python
new = pd.read_csv(f"{BUCKET}/.../summary_v4/{anc}.flexrv_vs_skato.v4.tsv.gz", sep="\t")
old = pd.read_csv(f"{BUCKET}/.../summary_v4_prefix/{anc}.flexrv_vs_skato.v4.tsv.gz", sep="\t")
j = new.merge(old, on=["anc", "trait", "gene"], suffixes=("_new", "_old"))
assert (j.source_new == "v7").any() and (j.source_old == "fixed").any(), "self-comparison"
b = j[(j.source_new == "v7") & (j.source_old.isin(["fixed", "v4"]))]   # BINARY genes
```

**`source_old == "fixed"` ALONE IS A BIASED SUBSET -- do not use it** (found
2026-10-01, after it had already produced a verdict). `flexrv_v4_fixed/` holds
only the chromosomes the EARLIER filtered rerun (#164/#166) touched; everywhere
else the summary falls back to `flexrv_v4/` and tags the row `source = v4`. So
that filter silently restricts the comparison to those chromosomes, and they were
selected by RECORDED SPA FAILURES -- a criterion ORTHOGONAL to #172's damage,
which is converged-but-wrong and was never recorded as a failure. Measured: it
cost eur 38 % of its genes (366,633 of 590,081, per-trait median 11,339 against a
max of 18,441) while leaving amr almost untouched, so it degrades one ancestry's
verdict and not another's. Both directories are pre-#172 states and both are valid
old sides; accepting either recovers full coverage. `source_new == "v7"` is what
identifies a binary trait, so quantitative rows stay out without naming them.

**Check the coverage before reading any verdict**, because the deficit is
invisible in the ratio:

```python
g = b.groupby("trait").size()
print(len(b), "genes | per-trait min", g.min(), "median", int(g.median()), "max", g.max())
```

A per-trait minimum far below the median means some traits are short. The NEW
side should be flat (the #172 rerun was `REGIONS=ALL`): eur reads 18,435-18,441.
A short NEW side would be a lost chromosome instead (LEDGER #170, read
`chrStatus`); a short OLD side is this.

Then: the identical fraction, the count moved more than 2x in each direction,
the median ratio, and the ten largest INCREASES with their `set_best`. On eas
93.8 % were identical, 92 genes moved more than 2x, and every one of the ten
largest increases was FemInf or MatHem out of fifteen traits -- chance would put
about 1.3 there -- which is how the seven spurious genes were found. **Rank by
the ratio, not by `p_skato / p_flexrv`**, and if you do use the latter exclude
`p_skato` above ~0.9 first: SKAT-O returns exactly 1 once its smallest per-rho p
exceeds ~0.78 (LEDGER #86), so the most NULL cells sort to the top as the most
suspicious.

**What a clean result means.** Identical everywhere is the expected outcome for
an ancestry with no tail models; it says the fix was surgical there, NOT that
the arm was never at risk. The damage needs the COMBINATION of a tiny fitted
probability and a weight one atom carries, so an ancestry whose models have no
tail has nothing for the transforms to concentrate on.


### 7b. The FlexRV-vs-SKAT-O discordance screen (2026-10-01)

Built and validated against a known positive; full result in
`docs/state/handoff-2026-10-01-flexrv-skato-discordance-screen.md`. Reads the two
gene-level summaries only -- no per-trait file, no new compute.

**The comparator is SKAT-O's `burden:` p at `@0.001`, not `skato:`, and not the
summary's `p_skato`.** FlexRV runs `--r.corr=1`, so it is a burden test at ONE max
MAF; `p_skato` pools all 18 cells including 0.01 and 0.0001, and `skato:` is a
variance-component test. Either mismatch manufactures discordance legitimately.

**The CCT must carry its tail branch.** A naive `0.5 - atan(S)/pi` floors at
`2^-54 = 5.551115e-17` and fabricates discordance wherever the comparator is
genuinely significant -- which, because quantitative traits have far more such
genes, produced a fake 20-vs-0 binary/quantitative asymmetry that read like a new
artefact. **If `5.551115e-17` appears in any CCT output in this repo, that is the
bug and not a p-value.**

```python
EPS = 1e-320
def cct(row):
    p = np.asarray(row, float); p = p[np.isfinite(p)]
    if p.size == 0: return np.nan
    p = np.clip(p, EPS, 1 - 1e-16)
    t = np.where(p < 1e-15, 1.0/(p*np.pi), np.tan((0.5 - p)*np.pi))
    S = t.mean()
    return 1.0/(np.pi*S) if S > 1e15 else 0.5 - np.arctan(S)/np.pi

MASKS = {"tight3": ["pLoF","damaging_missense","pLoF;damaging_missense"],
         "tight4": ["pLoF","damaging_missense","other_missense","pLoF;damaging_missense"]}
# cols are burden:<mask>@0.001 in <anc>.skato_cells.v4.tsv.gz; SAIGE writes
# combined masks with ';' although the command line takes ':'. ASSERT the column
# match: a wrong separator gives an empty list, NaN everywhere, and reads as
# "no discordance" rather than as a bug.
```

Then `orders = log10(p_burden_cauchy) - log10(p_flexrv)`, with **`p_burden > 0.9`
excluded** (LEDGER #86: SKAT-O returns exactly 1 once its smallest per-rho p
exceeds ~0.78, so the most NULL genes otherwise sort to the top as the most
suspicious). Report tight3, tight4 and the `min` over the same cells together; if
the flagged set moves between them the screen is too fragile to act on (found in
exactly one place: eur VTE `ENSG00000115718`, 4.89 on tight4 against 11.25 on
tight3).

**Two controls, and run them first.**
- POSITIVE: `amr EFRMB ENSG00000213030` (`p_flexrv` 2.81e-16, above the
  13-carrier cap, LEDGER #176) MUST appear at the top of amr. It reads 13.91 /
  13.78 / 13.34 on tight4 / tight3 / min. If it does not flag, fix the screen
  before believing anything else it says.
- NEGATIVE: the quantitative arm, which cannot carry the #172 artefact
  (`traitType == "binary"` gates the routing and the saddlepoint is binary-only).
  **CONFOUNDED BY POWER** -- quantitative traits are far better powered, so
  genuine FlexRV wins are commoner, and on this data the confound was the whole
  story. Do not read a binary-minus-quantitative difference as the artefact's
  footprint without checking what the quantitative hits ARE.

**What a real signal looks like, so a future run can tell them apart.** eur's
twenty most discordant quantitative genes are the canonical ones for their traits
-- GPT for ALT, GOT1 for AST, LDLR/APOB/PCSK9 for LDLC and TChol,
ABCA1/LCAT/APOA1 for HDLC, APOC3/ANGPTL3 for TG -- with `min` above 10 orders for
all twenty, only 4 of 20 on a pseudo-marker-only set, across seven traits. A
genuine FlexRV win is spread, survives the `min`, and names a gene that belongs
to the trait. An artefact concentrates, and LEDGER #176's case sits at 13.9
orders in a single gene of 1.65 million.


### 7c. Per-trait calibration as a COLUMN (2026-10-01)

`harness/tools/aou_flexrv_calibration.py` -- the third of the three
detectability instruments `docs/state/next.md` calls blocking. The null-median
check used by hand on 2026-09-30 is three lines, which is why it was a notebook
cell and why it was only ever run on the traits someone thought to look at.

```sh
python3 aou_flexrv_calibration.py {anc}.flexrv_vs_skato.v4.tsv.gz ... > calib.tsv
```

One row per (ancestry, trait, ARM) -- both `p_flexrv` and `p_skato` -- with
`n_genes`, `min_p`, the gene and set/cell it came from, and:

| column | meaning |
|---|---|
| `null_median` | `1 - 0.5**(1/m)`, the median of the smallest of m uniforms |
| `ratio` | `min_p / null_median`; 1.0 is exactly on the null, <1 extreme, >1 conservative |
| `p_cal` | `P(min <= obs) = 1 - (1-min_p)**m`, the calibration p-value |

**Why the MINIMUM and not a QQ slope or a lambda.** Every FlexRV defect found so
far (#163/#164, #172, #176) produces ONE absurd gene while leaving the bulk
untouched -- which a slope or a genomic-control factor averages away. The
minimum is the statistic that sees exactly that shape. (Separately, lambda is
invalid here anyway: burden p-values have atoms.)

**The discriminator, and the reason both arms are computed.** A real association
also gives a small minimum. One arm extreme with the other null is the shape all
three known defects had; BOTH arms extreme on the same gene is probably biology.
The table makes that readable without recomputation.

**A flagged trait is not wrong, it is worth looking at.** Apply multiplicity to
`p_cal` across traits before acting -- at a thousand traits, 0.1 %% of them below
`p_cal` 0.001 is the expectation, not a finding.

Self-test `aou_flexrv_calibration_selftest.sh` on fabricated summaries: the
closed forms at a RELATIVE tolerance (the table is 6 significant figures, and
the first version compared at 1e-12 absolute and failed for that reason alone),
20 NULL traits x 2 arms as the negative control, a planted FlexRV-only gene that
must fire on one arm and stay quiet on the other, and a planted concordant gene
that must fire on both.


## 8. Converged-but-wrong saddlepoints: `SPA_DUMP_BELOW` (2026-09-28, LEDGER #172)

The scan (section 1) flags FAILED saddlepoints (`Is.SPA FALSE`); a saddlepoint
that converges to a wrong value passes it, and the splice keeps its v4 row.
Three FlexRV cells did exactly that: amr BenCervUterNeo ENSG00000101198
`AM_sbeta20_MACsparse` p 2.2e-53 (normal p 8e-101), amr BenIntNeo ENSG00000196109
(sbeta sets), eas MatHem ENSG00000128607 `AM_sbeta10_MACsparse` p 3.5e-93 with a
NORMAL p of 0.02. Same carriers, sane p under `sbeta5` and `lof`; the trigger
is the weight spread of `sbeta10`/`sbeta20`. Every kernel change since the v4
image (`f2fb24da`) is on a failure path, so the current image is predicted to
reproduce the values -- the run has to be watched, not trusted.

With `SAIGE_SPA_DUMP_DIR` set, `SAIGE_SPA_DUMP_BELOW=<p>` ALSO dumps every
CONVERGED saddlepoint whose p is below `<p>`, as `spaBelow_<pid>_<k>.tsv`.
Every dump (both kinds) now carries a second header line -- `kernel`
(`fast`/`dense`), both roots and iteration counts, the q-side saddle's
validity, the Chernoff bound, `NAmu`/`NAsigma`, `var2`, the variance ratio,
the reported `p` and `pNA` -- and an `nz` column flagging the carriers, so
`aou_spa_dumpcheck.sh` replays the run's OWN kernel on the file: the `RUN:`
line is what the run saw, `FAST replay` prints `t*q-K`, `w`, `v`, `Ztest` at
the root (for these vectors the Legendre term is ~6 and `w` ~3.5; a `w` near
15 says K is ~115 too small at the root, or q is not the vector's), and
`MONTE CARLO reference` (2e6 draws, carriers exact, the rest Gaussian) is the
p the cell should have printed. Dumps from an older image (no second line,
no `nz`) still parse.

Image with the switch: `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-f14d47e5` (IMG8; the
same bytes on Docker Hub as `astheeggeggs/saige-slim:1.5.2-dev-f14d47e5`). Wrapper: inside the `SPA_DUMP` block,
`[[ -n ${SPA_DUMP_BELOW:-} ]] && export SAIGE_SPA_DUMP_BELOW="${SPA_DUMP_BELOW}"`.
Notebook, the three (trait, gene) pairs only, explicit regions rather than
the scan's:

```python
IMG8 = "gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-f14d47e5"
def submit_below(anc, pheno, gene, image=IMG8, thr="1e-10", trait="binary", **kw):
    return submit(name=f"fxbelow-{anc}-{tag(pheno)}", image=image,
                  masks="pLoF:damaging_missense:other_missense", fasttest="auto",
                  outdir=f"{anc}/flexrv_below", anc=anc, pheno=pheno, trait=trait,
                  chrs=range(1, 23), plk=lambda c, a=anc: cod(c, a), grp=am_grp,
                  pooled=True, sidecars=True,
                  env={"FLEXRV_SCORE": "AM", "FLEXRV_MAXMAF": "0.001", "MAXMAF": "0.001",
                       "RCORR": "1", "REGIONS": gene, "SPA_DUMP": "1", "SPA_DUMP_BELOW": thr}, **kw)
# masks/MAXMAF/model/cutoff exactly the v4 job's (RUNS["fxv4-..."]); then
# aou_spa_dumpcheck.sh on each .spaDump.tgz as in section 5
```

The threshold `1e-10` keeps the dump to the cells of interest (a real 1e-10
in one gene is at most a handful of sets). NOT to do yet: rerun every cell
that triggered the SPA. Nothing changes until the kernel bug is found; the
affected population is then binary traits only, and the lower bound from
output, egress-safe as a count per ancestry x trait, is rows with `Is.SPA
true` and `log10(p.value.NA) - log10(p.value) > 3`; the BenCervUterNeo
shape (p between the normal p and the truth) is not screenable, so the
measurement after the fix is a full binary FlexRV rerun diffed against v6.

## 9. LEDGER #172 FIXED: the AoU rerun of the binary FlexRV arm (2026-09-29)

Both mechanisms are on the row; the fix routes a custom or transform set whose
weight one atom carries to the exact test (Kish effective carrier count at or
below `--max_MAC_for_ER`, at most 13 material carriers), and every routed cell
is named in `<prefix>.exactByWeight.txt`. Which v6 cells were wrong cannot be
read off the output (mechanism A leaves a p between the normal p and the truth),
so the population is measured by rerunning and diffing, not screened:

1. **Image**: IMG9 = `gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-a08fb295` (the fix
   commit; the same bytes on Docker Hub as `astheeggeggs/saige-slim:1.5.2-dev-a08fb295`).
2. **Rerun every binary trait's FlexRV job of every ancestry, unfiltered** (no
   `REGIONS`): the atom-dominated cells are wherever a gene has one high-score
   rare variant, which is everywhere. Quantitative traits never touch the
   saddlepoint and are not rerun. Output directory `flexrv_v7/`. The wrapper
   needs no change; the new sidecar is delivered like `.spaFallbacks.txt.gz`
   if the wrapper globs `<prefix>.*` for its deliveries -- check one job's
   listing for `.exactByWeight.txt.gz` before launching the rest.
3. **Stage check** as in section 6: `SPA fallbacks: 0` everywhere, no lost
   chromosome, and now `Exact test by weight: N cell(s)` per job; the sum over
   jobs of N is the affected population, egress-safe as a count.
4. **Splice** over `flexrv_v4_fixed/` into `flexrv_v4_fixed2/` (or replace):
   `aou_spa_splice.sh` accepts a changed cell when the rerun's
   `.exactByWeight.txt.gz` lists it (Region|Weight), carries the sidecar
   restricted to the spliced regions, and prints
   `cells routed to the exact test by weight (LEDGER #172): N`. With an
   UNFILTERED rerun, `REGIONS` is every region of the trait (`cut -f1` of the
   v4 file minus `Cauchy`), and the control's other rules apply as before. Any
   refusal is read, not re-submitted.
5. **Controls**: (a) a routed cell's new p must sit in the decade of its own
   `lof`-only or `unw` twin when that twin took the exact test by MAC -- the
   set IS that atom; (b) the three cells that opened the row: BenCervUterNeo
   and BenIntNeo (amr), MatHem x ENSG00000128607 (eas) -- their sbeta10/sbeta20
   rows must now read ~1e-3 (amr) and ~1e-1 (eas), not 1e-53/1e-47/1e-93;
   (c) `Exact test by weight: 0` on every quantitative job if any is rerun by
   mistake (the routing is binary-only).
6. Then the summary and QQ jobs of section 7 over the spliced files.

**What the first sas pass found (2026-09-29).** Six jobs, `STAGE OK`, 0
fallbacks, 22 chromosomes each; ~570k-620k routed cells per trait (3.49 M over
the six, a third of all cells, as on the one chromosome inspected); and **no
`.exactByWeight.txt.gz` delivered**: the wrapper gathers only the sidecars it
names. Two edits fix that, then the jobs are resubmitted into the same
`flexrv_v7/`:

- wrapper `saige_step2_brava_main_original.sh`, the sidecar loop gains a pair:
  `"exactByWeight:${OUT_FILE_EXACTBYWEIGHT:-}"` beside
  `"spaFallbacks:${OUT_FILE_SPAFALLBACKS:-}"` (the loop body already gathers the
  per-chromosome parts or writes a placeholder);
- `submit()`, the sidecars block gains
  `"--output", f"OUT_FILE_EXACTBYWEIGHT={od}/{pre}.exactByWeight.txt.gz",`.

**LEDGER #176's counter needs the SAME two lines and does not have them
(2026-10-02).** `submit()` declares no `OUT_FILE_EXACTBYWEIGHTABOVECAP` and the
wrapper's deliver loop has no pair for it -- read in the source, so the wiring is
certainly absent and will bite the moment a count IS non-zero.

**What is NOT established: that anything was lost on the first job.** The sidecar
did not appear in that run's outputs, and the file is written ONLY when the count
is non-zero, so its absence is consistent with EITHER the missing wiring OR a
count of ZERO -- and a zero would be a FINDING about #176, not a plumbing fault.
The counter fires only when BOTH `n_eff <= --max_MAC_for_ER` AND
`n_material > 13`; a cell with 21 material carriers but `n_eff` of 8 is not
concentrated by that test and is correctly not counted, which is exactly the
"two different tests, so different p-values are right" outcome. **The job log
distinguishes them in one line and nothing else does.** Wire it anyway:

- `submit()`:
  `"--output", f"OUT_FILE_EXACTBYWEIGHTABOVECAP={od}/{pre}.exactByWeightAboveCap.txt.gz",`
- the wrapper: `OUT_FILE_EXACTBYWEIGHTABOVECAP="${OUT_FILE_EXACTBYWEIGHTABOVECAP:-}"`
  in the defaults, and `"exactByWeightAboveCap:${OUT_FILE_EXACTBYWEIGHTABOVECAP:-}"`
  in the `for pair in ...` deliver loop.

The file is written ONLY when the count is non-zero, so the loop must tolerate
its absence the way `.pooledTests` already does. Until that lands, the COUNT is
still in the job log (`Concentrated cells above the exact test's enumeration
cap: N`) and the per-cell `n_eff`/`n_material` can be had instead from a rerun
with `SPA_DUMP=1` and `SPA_DUMP_BELOW`, computing Kish's count from the dumped
vector in the workbench -- the dump is individual-level and does not leave.

**The control for a COMPLETE rerun** is the splice tool in its all-regions
mode: `REGIONS=ALL` (every region of the new results file, derived inside the
job) and `FAILURES=NONE` (the old side already carries the earlier fixes). It
refuses any changed cell that is neither listed in the rerun's
`.exactByWeight.txt.gz`, nor a MAC <= 4 cell (v4 rows never rerun took the
saddlepoint there before #166), nor an unobservable burden cell, and writes the
verified files to `flexrv_v7_checked/` under the v4 names:

```python
def submit_check_v7(v4name, image="gcr.io/google.com/cloudsdktool/google-cloud-cli:slim"):
    r = RUNS[v4name]; anc = r["anc"]; pre = r["out_prefix"].rsplit("/", 1)[-1]
    fixed = f"{BUCKET}/SAIGE_SLIM_TESTING/{anc}/flexrv_v4_fixed/{pre}"
    has_fixed = subprocess.run(["bash", "-lc", f"gcloud storage ls {fixed}.txt.gz"], capture_output=True).returncode == 0
    old = fixed if has_fixed else r["out_prefix"]          # the current final: spliced if it was, else v4
    new = f"{BUCKET}/SAIGE_SLIM_TESTING/{anc}/flexrv_v7/{pre.replace('fxv4', 'fxv7')}"
    args = ["aou_dsub", "--name", v4name.replace("fxv4", "check7"), "--image", image,
            "--env", f"OLD_PREFIX={old}", "--env", f"NEW_PREFIX={new}", "--env", "REGIONS=ALL", "--env", "FAILURES=NONE",
            "--output-recursive", f"OUT_DIR={BUCKET}/SAIGE_SLIM_TESTING/{anc}/flexrv_v7_checked",
            "--machine-type", "n2-standard-2", "--disk-size", "20", "--script", "aou_spa_splice.sh"]
    script = "source /home/jupyter/aou_dsub.bash\n" + " ".join(shlex.quote(a) for a in args)
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=os.environ)
    print(v4name, "old =", "fixed" if has_fixed else "v4", "|", (out.stdout or out.stderr).strip().splitlines()[-1])
```

The log's lines to read: `cells routed to the exact test by weight (LEDGER
#172): N`, the unobservable-cell line with its direction counts, and the
absence of `ERROR`, and -- since 2026-09-30 -- `rows accepted on statistic
identity alone (CONTROL 1c, the weaker tier): N` in `<pre>.splice.txt`, which
is NOT a pass but a count of rows whose cell the control could not name (14b).
The notebook's `aou_spa_splice.sh` must be the repo's current version; its
self-test prints what it ran, so compare that line rather than a number
remembered from here (it was 3 positive controls and 7 planted refusals before
2026-09-30, and 8 and 9 after -- a stale count in this file is exactly the
trap section 14c describes).

`SAIGE_EXACT_BY_WEIGHT=0` in a job's env restores the pre-#172 routing for an
A/B on one trait, which is how to measure the population on a sample before
committing to all of it.

## 10. How many samples are pinned, per model (LEDGER #174, 2026-09-29)

Every step-2 run since the fit gates prints, when any fitted value is below
1e-8 or above 1 - 1e-8, `the null model in <file> has K of N fitted values (x%)
at the 0/1 boundary (range a .. b)`. `harness/tools/aou_pinned_scan.sh` greps
the newest delivered log of every job with a given name prefix and writes ONE
table, one row per (ancestry, trait), with any count from 1 to 19 written as
`<20` and its percentage masked. No model is opened and nothing is rerun.
Self-test: `aou_pinned_scan_selftest.sh`.

```python
def submit_pinned(prefix="fxv4-", tag="v4", image="gcr.io/google.com/cloudsdktool/google-cloud-cli:slim"):
    args = ["aou_dsub", "--name", f"pinned-{tag}", "--image", image,
            "--env", f"LOGS={BUCKET}/dsub/logs", "--env", f"PREFIX={prefix}", "--env", f"TAG={tag}",
            "--output-recursive", f"OUT_DIR={BUCKET}/SAIGE_SLIM_TESTING/pinned",
            "--machine-type", "n2-standard-2", "--disk-size", "20", "--script", "aou_pinned_scan.sh"]
    script = "source /home/jupyter/aou_dsub.bash\n" + " ".join(shlex.quote(a) for a in args)
    r = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=os.environ)
    print((r.stdout or r.stderr).strip().splitlines()[-1])
```

Controls inside the table: `eas`/`mathem` must read 77 of 4,089 with a smallest
fitted value of 2.22e-16 (measured from its dump job's log on 2026-09-29), and
the six sas binary traits must read `none_reported` with gates `on_or_unknown`
(measured the same day from the v7 logs). If the eas row reads `none_reported`,
the logs of that prefix predate the gates or ran with them off: use the prefix
of a later run (`fxv5-`, `fxv7-`).

## 11. Where do the models' fitted values sit? (LEDGER #174/#175, 2026-09-29)

The pinned-sample guard's threshold was set at 1e-8 on the assumption that a
separated group is fitted below it. The bundled fixture says it is not: the
group is left wherever glm stopped, and straddles the threshold. The step-2
logs report only the count at or below 1e-8 (section 10), so this reads the
models themselves.

`aou_mu_decades.sh` runs on the SAIGE image (it needs R and nothing else).
dsub localizes the models inside the workspace, so nothing is egressed. The
table has one row per binary model: the number of samples in each decade of
flip probability, the count and the mass at or below 1e-8, 1e-7, 1e-6 and
1e-5, and the widest run of empty decades. Every count from 1 to 19 is
written `<20` and the mass that goes with it `masked`. It is safe to paste.

Upload `aou_mu_decades.sh`. The models of four ancestries share ONE directory,
`data/saige/run/nullglmm_v4/<anc>_<TRAIT>.rda`, and amr's are in
`nullglmm_v4_rc0.125/` (`docs/state/aou-saige-parameters.md`); they are NOT
under `SAIGE_SLIM_TESTING/`, which is where the first version of this recipe
looked and found nothing. So the job takes a file pattern, and dsub localizes
that ancestry's models only. Locate them first; this prints counts and nothing
else:

```python
MODEL_DIR = {a: f"{BUCKET}/data/saige/run/nullglmm_v4" for a in ["afr", "eas", "eur", "sas"]}
MODEL_DIR["amr"] = f"{BUCKET}/data/saige/run/nullglmm_v4_rc0.125"
try:
    print("the notebook's STEP1:", str(STEP1).replace(BUCKET, "<bucket>"))
except NameError:
    print("STEP1 is not defined in this session")
for anc in ["eas", "sas", "afr", "eur", "amr"]:
    n = subprocess.run(["bash", "-lc", f"gcloud storage ls '{MODEL_DIR[anc]}/{anc}_*.rda' 2>/dev/null | wc -l"],
                       capture_output=True, text=True).stdout.strip()
    print(anc, "model files:", n, "(expected 25 / 15 / 42 / 43 / 42 for eas / sas / afr / eur / amr)")
```

Confirmed 2026-09-29: `STEP1` names exactly these two directories and the
counts read 25 / 15 / 42 / 43 / 42. If a count is ever 0 and `STEP1` names
another directory, put that directory in `MODEL_DIR`. Then:

```python
def submit_mu_decades(anc, image=IMG9):
    args = ["aou_dsub", "--name", f"mudec-{anc}", "--image", image,
            "--input", f"MODELS={MODEL_DIR[anc]}/{anc}_*.rda", "--env", f"TAG={anc}",
            "--output-recursive", f"OUT_DIR={BUCKET}/SAIGE_SLIM_TESTING/mu_decades",
            "--machine-type", "n2-standard-4", "--disk-size", "100", "--script", "aou_mu_decades.sh"]
    script = "source /home/jupyter/aou_dsub.bash\n" + " ".join(shlex.quote(a) for a in args)
    r = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=os.environ)
    print((r.stdout or r.stderr).strip().splitlines()[-1])
submit_mu_decades("eas")
```

The job fails, with `ERROR: no model file matched` in its log, if the pattern
matches nothing: a wrong path cannot look like a clean result.

Start with eas: it is the smallest and holds two of the seven models. Reading
the result:

```python
t = load(f"{BUCKET}/SAIGE_SLIM_TESTING/mu_decades/mu_decades.eas.tsv.gz")
cols = ["model", "N", "n_le_1e-08", "n_le_1e-07", "n_le_1e-06", "n_le_1e-05",
        "mass_le_1e-08", "mass_le_1e-06", "empty_decades_in_widest_gap"]
print(t[cols].to_string(index=False))
```

What to read:

- **`n_le_1e-06` against `n_le_1e-08`.** Equal: the model's separated samples
  are all below the shipped threshold and the guard is complete for it. Larger:
  the group straddles the threshold, as on the fixture, and the guard is partial.
- **A model with `n_le_1e-08` of 0 and `n_le_1e-06` above 0.** A separated group
  that the log scan could not see (LEDGER #175).
- **`empty_decades_in_widest_gap`.** Three or more empty decades between the
  lowest group and the rest is what separation looks like; a smooth tail down
  to 1e-6 is a model that really does assign tiny probabilities, and those
  samples are not separated.
- **`mass_le_1e-06`.** The most any one test could leave out at that threshold.

## 12. Which covariate carries the tail? (LEDGER #174/#175, 2026-09-29)

The eas decade scan (section 11) found a tail in two of fifteen binary models,
FemInf and MatHem, and it is NOT the shape of the fixture's separated group:
the widest gap is 0 decades and the counts rise steadily, 62 / 103 / 175 / 310
and 77 / 141 / 244 / 394 at or below 1e-8 / 1e-7 / 1e-6 / 1e-5. A categorical
covariate that separates a group gives a cluster with a gap above it; a
continuum is what a steep CONTINUOUS effect gives. Both traits have N near
4,050, like the three other female-specific traits (3,961 to 4,098) and unlike
the both-sex traits (6,090 to 6,379), so the other sex was excluded.

`aou_mu_drivers.sh` answers the rest from the model alone: the covariates that
entered the fit (a constant or collinear column is dropped by step 1 before
the fit and is simply absent), the case count, and which covariate carries the
tail samples' linear predictor away from everyone else's. Every number is a
sum over the cohort or over the tail group; a tail of 1 to 19 samples is
written `<20` and its rows are withheld. Self-test: 16 checks on simulated
models where the driver is known (a separating covariate; age and its square).

Upload `aou_mu_drivers.sh`, then (BreastCanc is the control: female-specific,
no tail, so its summary is withheld but its covariates and case count print):

```python
def submit_mu_drivers(anc, traits, image=IMG9):
    for tr in traits:
        args = ["aou_dsub", "--name", f"mudrv-{anc}-{tr.lower()}", "--image", image,
                "--input", f"MODELS={MODEL_DIR[anc]}/{anc}_{tr}.rda", "--env", f"TAG={anc}_{tr}",
                "--output-recursive", f"OUT_DIR={BUCKET}/SAIGE_SLIM_TESTING/mu_drivers",
                "--machine-type", "n2-standard-4", "--disk-size", "50", "--script", "aou_mu_drivers.sh"]
        script = "source /home/jupyter/aou_dsub.bash\n" + " ".join(shlex.quote(a) for a in args)
        r = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=os.environ)
        print(tr, (r.stdout or r.stderr).strip().splitlines()[-1])
submit_mu_drivers("eas", ["FemInf", "MatHem", "BreastCanc"])
```

Reading it, from the job logs:

```python
for tr in ["FemInf", "MatHem", "BreastCanc"]:
    log = fetch_log(f"mudrv-eas-{tr.lower()}")
    print(log[log.find("====="):][:1600])
```

What each answer looks like:

- **`columns whose name contains 'sex': none`** on all three: sex did not enter
  these fits (not passed, or passed and dropped as constant).
- **a sex column present**: it survived the collinearity check, so the cohort
  is not single-sex after all. Then look at how much of the gap it carries.
- **`age` and `age2` carry the gap**: a steep age effect. The samples in the
  tail are those the model places furthest along it.
- **one categorical column carries the gap**: a separated level, the fixture's
  shape, whatever the decade table's gap column said.
- **`(no covariate: random effect)` carries the gap**: the tail is the random
  effect's doing, and the covariates are not the story.

## 13. The eas answer, and the four checks it sets up (2026-09-29)

Section 12's three jobs came back. `age` and `age2` carry **99.8 %** of the
tail's gap on eas FemInf (63.3 + 36.5) and **100.0 %** on eas MatHem
(50.1 + 49.9); no PC exceeds 0.2 %, the random effect contributes -0.1 %, no
categorical column appears, and `sex` entered neither fit. BreastCanc -- the
control: same twelve covariates, same female-only cohort, N 4,074 against 4,054
and 4,089 -- has a tail of **0**. So the `[ASSUMPTION: the tail is a steep age
effect. Not measured.]` of `handoff-2026-09-29-step-back.md` is MEASURED, and
the covariate design, the PC set and the sex restriction are all ruled out as
the cause: what differs is the trait's relationship with age. FemInf and MatHem
can only occur in reproductive-age women and the model puts everyone else at
~0; breast-cancer risk rises monotonically, so nobody is put at ~0.

It reconciles two ways, which is the check to make before believing it.
FemInf's case fraction is 103/4,054 = 2.54 %, logit -3.65; minus the 13.12 gap
that is `exp(-16.8)` = 5.6e-8. MatHem: -3.66 - 13.0 -> 5.8e-8. Both land
between the 1e-8 and 1e-6 rungs, which is where section 11's table put 62/175
and 77/244 samples.

**This IS separation, on a CONTINUOUS covariate** -- what `R/FitGates.R:249`
says the mu check exists to catch. What is absent is the CATEGORICAL shape, so
the covariate-cell scan is not the relevant half and the fixture's
cluster-with-a-gap is not what real data does here.

**The census that produced "the seven" is at the wrong threshold.** It is
section 10's pinned scan, whose count is the gate's own `eps = 1e-8`
(`checkMuSeparation`), and the measured tails sit ABOVE that: only 62 of
FemInf's 175 tail samples are below 1e-8. A model with a slightly shallower age
effect keeps its whole tail between 1e-8 and 1e-6 and never appears in the
seven. Section 11's decade scan, which counts at 1e-5 too, has been run for
**eas only**.

### 13a. Recount, the other four ancestries (do this first)

One `mudec-<anc>` job per ancestry, model-only, no genotype pass, nothing
egressed. `submit_mu_decades` is section 11's.

```python
for anc in ["sas", "afr", "eur", "amr"]:
    submit_mu_decades(anc)
```

Flag every model whose count at or below **1e-6** is not `<20`, not just those
with something at 1e-8. A trait that is flagged at 1e-6 and absent from the
seven is the finding this section exists for.

### 13b. Pre-registered prediction (written 2026-09-29, BEFORE 13a was run)

The susceptibility rule the eas answer implies: a trait is susceptible when its
probability falls to ~0 over part of the cohort's OWN age range, steeply enough
to reach 1e-6. "Sex-specific" is not enough and neither is "age-dependent" --
BreastCanc is both and has no tail, because its risk is monotone rather than
confined to a window.

- FemInf and MatHem: age-driven tail in EVERY ancestry they were fitted in.
- amr AMD: age-driven with the **sign reversed** -- its tail should be the
  YOUNG, since AMD risk climbs with age.
  `[ASSUMPTION: AMD is age-related macular degeneration, read from the
  abbreviation, not verified against the phenotype definition.]`
- Any other pregnancy- or childbirth-confined trait in the list: tail.
- Monotone sex-specific traits of BreastCanc's shape: NO tail.

A flagged trait that is on none of those lines falsifies the rule, and that is
the outcome worth having.

### 13c. The gate's own record -- ALREADY MEASURED, one number still to read

No new job: section 10's warning IS `checkMuSeparation`'s, so the gate already
recorded this at fit time and the pinned table holds it. Its control line
reads eas/mathem at **77 of 4,089** with a smallest fitted value of 2.22e-16 --
1.88 %, below `--separationMaxFraction=0.02`, hence a WARNING and not a
refusal. **Prediction for the row not yet quoted: eas FemInf reads 62 of 4,054
= 1.53 %, also a warning.** If it reads anything else, the arithmetic above is
wrong and 13a's reading needs redoing. Both models passed by a whisker; counted
at 1e-6 they are 4.3 % and 6.0 % and both would be refused -- which is LEDGER
#175 exercised by real data, not by the fixture.

A stricter threshold would NOT have been the fix: `FitGates.R` refuses the
whole fit, so there is no per-gene or per-sample granularity -- it would have
cost each affected trait entirely, all 22 chromosomes, and a `stop()` takes a
multi-phenotype batch with it.

### 13d. Which weight sets drive an extreme gene

Weight concentration is the standing explanation for why the SKAT-O arm is
quiet while FlexRV is not: the hazard needs a tiny-mu carrier AND enough weight
on that one carrier for the root to pass `-log(mu)`, and beta(1,25) mask
weights never concentrate that far (which is why the beta-weighted `region`
TSAN workload never enters the #172 routing). It predicts WHICH rows are
extreme, and the per-set rows are already in the FlexRV output -- no rerun.

```python
d = load(f"{BUCKET}/SAIGE_SLIM_TESTING/<anc>/flexrv_v4/<file>")
g = d[d.Region == "<GENE>"].sort_values("Pvalue_Burden")
print(g[["Weight", "Pvalue_Burden", "BETA_Burden", "SE_Burden", "MAC",
         "Number_rare", "Number_ultra_rare"]].head(20).to_string(index=False))
```

**Prediction: the extreme rows are the CONCENTRATED sets** -- `sbeta20`, the
steep `cbrt` rungs, and the `MAC1`/`MAC3` sets on cells whose `Number_rare` is
0 (the collapsed pseudo-marker alone). If `unw_unw` and the gentle beta sets
are extreme too, weight concentration is not the whole story and the
discordance signature should not be leaned on for the other four ancestries.

### 13e. The damage check needs no special rerun

Whether the seven actually produced bad results is section 9's gene-level
comparison on the IMG9 rerun, and amr/afr/eur are **not started**. Run them for
ALL traits, not a susceptible subset: the #172 fix moves results wherever a
one-atom set exists, so the rerun is needed regardless, and it delivers the
damage check for the seven as a by-product. A subset rerun would cost a second
pass later and would also make the comparison's denominator a chosen set rather
than the whole arm.

## 14. The scan stamp is PER ANCESTRY, and a wrong one reads as clean (2026-09-30)

Cost most of an hour on afr. The notebook's `STAMP` is a SCALAR holding one
ancestry's stamp; the scans were minted in a submission loop seconds apart, so
afr's is `20260926-170525` while `STAMP` held `20260926-170535` -- **ten seconds
and one digit apart**.

**Every check downstream of it then returns "clean" from an empty directory.**
`regions_for` runs `gcloud storage cat` and returns `[]` on failure, and `[]`
means "nothing to rerun", i.e. BENIGN. A `gcloud storage cat` of a non-existent
`failures.tsv` likewise gives 0 rows through `awk`. So a one-digit stamp error
reported eleven confident BENIGN verdicts and zero failure rows, and nothing in
the output looked wrong.

**Derive the stamp, never trust the variable**, and make the count a hard gate:

```python
def stamp_for(anc, arm="flexrv_v4"):
    p = f"{BUCKET}/SAIGE_SLIM_TESTING/spa_scan/{anc}/{arm}/"
    out = subprocess.run(["bash", "-lc", f"gcloud storage ls '{p}'"],
                         capture_output=True, text=True).stdout.split()
    st = sorted(x.rstrip("/").rsplit("/", 1)[-1] for x in out if x.startswith("gs://"))
    if len(st) != 1: print(f"  WARNING {anc}: {len(st)} stamps {st} -- using the latest")
    return st[-1] if st else None
```

then `assert` that `rerun/*.regions.txt` is non-empty and that `failures.tsv`
parses to a non-zero row count BEFORE classifying anything. The tell that the
path was wrong on afr was `regions files: 0 of 124` while 113 of 124 prefixes
HAD a `flexrv_v4_fixed/` copy -- copies that only a splice reading a regions list
could have produced.

**What the classification is for.** `submit_check_v7` hardcodes `FAILURES=NONE`,
which is only legitimate when the old side already carries the #163-#168 fixes.
A prefix WITHOUT a `flexrv_v4_fixed/` copy is benign only if it has BOTH no
regions file and no row in `failures.tsv`; either present means the earlier round
had work there and did not finish it, and that prefix needs its `failures.tsv`
passed instead of `NONE`.

**Why "no regions file" is trustworthy as a signal**: the rerun list is failures
UNION tail regions (hence `tail_regions.tsv` beside `failures.tsv`), so on afr
`failures.tsv` covers 95 prefixes while 155 have regions files. A prefix with no
regions file therefore had neither.

**Measured on afr (2026-09-30)**: 124 binary jobs over 31 traits in 4
size-balanced, INTERLEAVED chromosome batches (`b0` = 1, 6, 9, 13, 15; `b2` = 2,
3, 5, 8, 20, 21 -- so a min-max summary of a batch is misleading). 113 prefixes
have a fixed copy; the other 11 are all benign, so **afr takes `FAILURES=NONE`
like eas**.

### 14a. `submit_check_v7`'s `has_fixed` fallback is LOAD-BEARING off eas (2026-09-30)

Section 9's `submit_check_v7` chooses its old side as `flexrv_v4_fixed/<pre>` when
that file exists and `flexrv_v4/<pre>` otherwise. **On eas that branch never
fires** -- all 15 prefixes have a fixed copy -- so a simplified version that
hardcodes `flexrv_v4_fixed/` passes the entire eas chain and then dies on
localization for every prefix the earlier round found nothing to rerun in.

Measured: **130 of 444 check jobs failed that way** (122 of eur's 320, 6 of amr's
124; afr's 11 had not been submitted yet), every one of them a prefix with no
fixed copy and NONE failing for any other reason. The share is what makes this
bite off eas: a finely batched ancestry has many batches with nothing to rerun
(eur is split ~10 ways per trait against eas's 1), so 38 % of eur's prefixes take
the fallback.

`FAILURES=NONE` stays correct for exactly these prefixes -- they are the benign
ones, with no regions file and no row in `failures.tsv`, so the raw `flexrv_v4/`
file already carries the #163-#168 state.

The failure signature is a job dead in under 60 s with a missing-input message,
which is the dsub localization signature generally. Diagnose by joining the
failed job names back to `RUNS` and testing for the fixed copy: if every failure
lacks one, it is this and nothing else.

### 14b. The splice control's pseudo-marker twin rule has a gap (2026-09-30, FIRST refusal of the rollout; FIXED 2026-09-30, and the first proposed fix was WRONG)

`aou_spa_splice.sh` refused two amr prefixes -- `pancreat-b0` and `stroke-b0` --
with `140 changed row(s) are neither a recorded failure nor a MAC <= 4
cell/marker nor an exact-by-weight cell nor a pinned-sample test` on
`.singleAssoc.txt.gz` only; `.txt.gz` passed `all explained` on both. **This is
the control working, and the gap is in the control, not the data.**

Same gene in both traits (`ENSG00000226314`), same count (140), same sets --
`AM_cbrt0.2_MAC100`, `AM_cbrt0.2_MAC50`, `AM_cbrt0.2_mafbeta2`.

**The mechanism.** The refused rows are ultra-rare pseudo-markers whose change IS
explained by LEDGER #166. But #166 eligibility is `MAC_nw <= --max_MAC_for_ER`,
the set's UNWEIGHTED collapse count (`src/Main.cpp:4132-4134`), and **that number
is in no delivered file** when the cell has non-collapsed rare markers: the
results `MAC` column is `MAC_nw` PLUS the kept rare markers' MACs, and
`.markerList.txt` is per CELL, not per SET. So the rule infers it from the
"0/1-weight twin with the same zero pattern", here `AM_unw_MAC100`, **whose count
is 119**. The twin bounds the set's count from ABOVE, which is what the script's
comment correctly says -- but the rule then uses the twin as a proxy for "is this
set at or below the cutoff", which only works when the twin is itself small.

**Why the two counts differ so much, and it is NOT a score threshold.** From
`docs/design/flexrv-design.md:85`, `cbrt_t(s) = (sign(s-t)|s-t|^(1/3) +
t^(1/3)) / (|1-t|^(1/3) + t^(1/3))`, which is 0 at `s = 0` and positive for
every `s > 0`; line 342 records that 0 and 1 are fixed points of every
transform. So a nonlinear score transform drops only markers whose score is
EXACTLY 0 -- and in these group files a 0 means ABSENT FROM AlphaMissense, not
benign. This gene's ultra-rare markers are mostly unscored, so `unw` keeps all
119 alleles while every nonlinear transform keeps only the scored markers'
alleles. Expect this wherever a gene's rare variants are largely unscored.

#### The fix this section first proposed is WRONG. Do not implement it.

The original text said, quoted here on purpose:
<!-- docdrift-allow: the next sentence is the false claim this section refutes -->
`AC_Allele2` on a pseudo-marker IS the unweighted collapse (LEDGER #118), so the
rule should read the row's own folded `AC_Allele2` instead of inferring via the
twin -- one awk change. **Both halves of that are false**, checked 2026-09-30
against the source and against a committed golden.

A pseudo-marker's `AC_Allele2` is the **WEIGHTED** collapse:
`altCounts = cs.paritySumDense(genoURVec.memptr())` at `src/Main.cpp:3888`, where
`genoURVec` is `max_m(w_r(m) g_m)`. #118 left it weighted DELIBERATELY -- its own
STATUS says "the UR pseudo-marker's OWN single-variant test context still reads
the weighted MAC, because that row tests the weighted pseudo-marker; that is why
every `*.singleAssoc.txt` is byte-identical". The comment does not predate #118;
it survived it on purpose.

MEASURED on the committed golden `harness/goldens/step2/plink_group_flexrv_skato.txt`.
GENE1's four weight sets have **identical per-set unweighted composition** -- all
three columns are per-set since #118 -- yet their pseudo-marker rows disagree:

| Region\|Weight | MAC | Number_rare | Number_ultra_rare | pseudo `AC_Allele2` |
|---|---|---|---|---|
| `REVEL_unw_unw` | 183 | 6 | 17 | **54** |
| `REVEL_unw_mafbeta2@0.02` | 183 | 6 | 17 | **47.85** |
| `REVEL_sbeta2_unw` | 183 | 6 | 17 | **50.47** |
| `REVEL_sbeta2_mafbeta2@0.02` | 183 | 6 | 17 | **44.4398** |

Same markers, same alleles; `AC_Allele2` spans 44.44-54 purely on the weights
(GENE2 is worse: 24.517 against 41). The 0/1 set reports an integer, which is
*why* the twin rule works and why the row's own value cannot replace it.
Comparing that against an allele-count cutoff is **LEDGER #134's defect** -- the
SCALE of a group file's weights choosing the statistical treatment -- moved into
the control, and it is **ANTICONSERVATIVE**: every down-scaling transform
(`cbrt0.2`, `mafbeta2`, i.e. precisely the refused sets) reads low, so the
control would have accepted rows that are not #166-eligible, on every ancestry,
silently.

#### What landed instead: CONTROL 1c, statistic identity (2026-09-30)

A last tier, singleAssoc PSEUDO-MARKER rows only, reached only after every named
rule above has failed: **accept a changed row iff `SE` and `p.value` are the only
columns that moved.** Measured on `input/group_new_snpid_withWeights_atom.txt`,
no AoU data needed:

| arm | rows moved | columns that moved |
|---|---|---|
| POSITIVE, #166 (`--max_MAC_for_ER` 0 -> 50) | 8 | `SE`, `p.value` |
| POSITIVE, #172 (`SAIGE_EXACT_BY_WEIGHT` 0 -> 1) | 10 | `SE`, `p.value` |
| NEGATIVE, variance ratio perturbed 0.1% | 14 | `BETA`, `SE`, `p.value`, `p.value.NA`, `var` |

`BETA` and `Tstat` are the score statistic, `var` its variance, `AC_Allele2` the
genotype summary, and the exact test re-derives `SE` from the new p. **`var` is
the discriminating column, not `Tstat`** -- the score does not depend on the
variance ratio, so `Tstat` alone would have passed the negative control.

**The limit, and why A still follows.** This is a FINGERPRINT, not a REASON. It
says the statistic is untouched and only the p-value method moved; it does NOT
say the method chosen was right, so a cell that switched methods WRONGLY
(LEDGER #176's shape) would pass. Every accepted row is therefore counted on
stderr (`pseudo-marker rows accepted on STATISTIC IDENTITY ALONE`), written to
`<pre>.spliceStatId.txt` beside the spliced file, and totalled in
`<pre>.splice.txt`, so the weaker provenance travels with the data instead of
being invisible. A refusal now also prints `columns that differ old -> new`,
which is what turned this diagnosis from guesswork into one line.

#### PRE-REGISTERED: what the re-splice of the four blocked prefixes must report

Written 2026-09-30, BEFORE the four are re-spliced, so the result can falsify it.

CONTROL 1c is the LAST tier: it only ever sees rows that no named rule
explained. On amr `pancreat-b0` and `stroke-b0` that set was exactly **140
rows**. So each of those two prefixes must splice clean with

    rows accepted on statistic identity alone (CONTROL 1c, the weaker tier): 140

in `<pre>.splice.txt`, and 140 MarkerIDs in `<pre>.spliceStatId.txt`.

- **Exactly 140, clean splice** -- the gap was the control's, as diagnosed.
- **FEWER than 140, still refused** -- this is the outcome to hope for and to
  read carefully. It means some of those rows move `BETA`, `var` or
  `AC_Allele2`, which the exact test never touches, so they are a REAL
  discrepancy in that data and not a control gap at all. The refusal now prints
  `columns that differ old -> new` for the first five; read that line before
  proposing any cause.
- **MORE than 140** -- impossible unless the tier order is wrong (1c runs after
  every named rule). Treat it as a bug in the control, not a finding.

afr `Asth` and `VTE` have no recorded refusal count, so they are observation
rather than prediction; note whatever they report, since a count far above ~140
would say the gap is wider than the two amr prefixes suggested.

Whatever the four report, the honest statement for a row accepted here is **"no
statistic moved"**, NOT "the row is explained" -- see the limit stated above.

#### OUTCOME, 2026-10-01: all four spliced clean, and the prediction was 2x out

| prefix | rows accepted on statistic identity |
|---|---|
| `fxv4-amr-pancreat-b0` | **70** |
| `fxv4-amr-stroke-b0` | **70** |
| `fxv4-afr-asth-b0` | **84** |
| `fxv4-afr-vte-b0` | **84** |

All four wrote `.splice.txt`, which is only reached after every `splice_one`
succeeds -- so none refused. amr and afr are now 124/124.

**The prediction of 140 was right in substance and wrong by exactly a factor of
two, and the reason is in this script.** `changed.rows` is
`comm -3 o.all n.all`, which emits the lines unique to EACH side, so one modified
row appears TWICE -- its old version and its new version. The refusal message
counts LINES (`bad++` per line) and said 140; CONTROL 1c counts distinct
MarkerIDs (`sidrow[$3]`) and says 70. **When comparing a refusal count against an
acceptance count, halve the refusal count first.**

**Three structural signatures now support 14b's account**, none of which was
predicted and all of which would be odd coincidences otherwise:
1. All four refusals are the SAME batch, `b0` = `chr1_6_9_13_15`, in both
   ancestries.
2. The two amr traits give the SAME count (70) and the two afr traits give the
   SAME count (84) -- same genotypes and group file, different phenotype, so the
   affected rows are fixed by the ANNOTATION STRUCTURE and not by the trait.
3. The counts differ BETWEEN ancestries (70 vs 84), which is what differing
   variant sets after per-ancestry QC should give.

That is consistent with "a gene whose rare variants are largely unscored by
AlphaMissense, so `unw` keeps every allele while every nonlinear transform keeps
only the scored ones" -- the mechanism 14b names. It remains an account of WHY
the twin bound is loose, not proof that each row's new p is right; CONTROL 1c
says only that no statistic moved.

#### The exact fix, SPECIFIED and NOT built (owner's call, 2026-09-30)

CONTROL 1c is a backstop. The exact fix is for the RUN to disclose the #166
routing, as `.exactByWeight.txt` and `.spaPinned.txt` already do for #172 and
#174 -- that asymmetry is the whole reason this gap exists, since those two
explanations are read from a sidecar and only #166 is inferred.

**It must NOT be always-on, and that is a measurement, not a preference.**
On `harness/goldens/step2/plink_group_flexrv_burden.txt`, 32 of 224 cells
(**14.3 %**) are ultra-rare-only at MAC <= 4, i.e. exactly what such a sidecar
would list. Scaled to an AoU FlexRV run (~18,000 genes x 112-192 sets per gene)
that is millions of rows per trait x chromosome -- plausibly larger than the
results file it documents, on a pipeline that already watches its egress.

**Shape: env-gated, default OFF**, the precedent being `SAIGE_SPA_DUMP_DIR`:

- `SAIGE_EXACT_TEST_LOG=<path>`; unset writes nothing, so **no golden moves**
  and the gate is a proof of INERTNESS rather than a re-baseline. That is the
  cheapest correct shape available and is why it is the recommended one.
- Two call sites, both in `src/Main.cpp`, and both already compute the number:
  the ultra-rare pseudo-marker's own single test (`cm.isER = (MAC_nw <=
  g_cfg.MACCutoffforER && ...)`, ~`:4134`) and the BURDEN emitter's twin
  (~`:4585`). `recordExactByWeight` is the working model for the writer.
- Schema, mirroring `.exactByWeight.txt` so one reader serves both:
  `Kind  Region  Group  max_MAF  Weight  MarkerID  MAC_nw  Reason`, with
  `Reason` = `mac` here (`weight` is #172's file). `MarkerID` is what keys a
  pseudo-marker row in the singleAssoc, so the splice control can join directly.
- Turn it on for any rerun that will be SPLICED; that is the only situation
  needing it, and it is decided at submit time.

**Then in `aou_spa_splice.sh`**: read `<NEW_PREFIX>.exactTestCells.txt.gz` if it
exists, key `MarkerID`, and accept a listed row as a NAMED explanation ahead of
CONTROL 1c -- the same tier structure `ewcells.txt` and `pncells.txt` already
use. 1c stays as the fallback for files produced without the log, which is every
file delivered before this lands, including the four blocked prefixes.

**What this does NOT fix:** LEDGER #176. A cell that switched methods WRONGLY is
still explained by its own listing. #176 needs the counter for concentrated
cells above the 13-carrier enumeration cap, which is a different instrument and
is item (1) of `docs/state/next.md`'s blocking list.

Drilled in `aou_spa_splice_selftest.sh`: the `ur` and `twin` plants (which move
`p.value` alone) flipped from planted refusals to positive controls for the new
tier, and three new plants -- `sidvar`, `sidbeta`, `sidac` -- are the same rows
with `var`, `BETA` or `AC_Allele2` moved as well and must still be refused.
`sidvar` differs from the accepted `ur` in exactly one column, which is the A/B
that shows the tier discriminates rather than waving everything through.

**Do not `--force` past a refusal.** Nothing wrong was written; a refusal is the
control doing its job, and the four blocked prefixes were the only ones in the
whole rollout across sas, eas, afr, amr and eur.

### 14c. The `summary_v4/` tables were STALE, and sas's "quiet" verdict rested on one (2026-09-30)

`aou_flexrv_summary.sh` writes `summary_v4/<anc>.flexrv_vs_skato.v4.tsv.gz`
whatever it read -- the filename and `RUNTAG` are identical for a pre-fix and a
post-fix build -- so an existing file looks like a completed step and says
nothing about which inputs went into it.

Measured on 2026-09-30, object timestamps: afr, amr, eur and **sas** all dated
**2026-09-28**, i.e. BEFORE the #172 reruns; only eas was rebuilt (that day,
11:35). sas's table reads `{'v4': 163347, 'fixed': 109589}` -- **not one `v7`
row**.

**So the step-back handoff's record of sas as "done, quiet" for summary and
comparison was drawn from PRE-FIX data.** sas's stage check and splice control
are independent of the summary and stand; it is only the GENE-LEVEL verdict that
is unverified, and `submit_summary("sas")` must be re-run before it can be
claimed. The same applies to afr, amr and eur: their 09-28 tables are pre-fix and
must be rebuilt AFTER their splices complete, not before -- a rebuild over an
incomplete `flexrv_v7_checked/` silently falls back to the pre-fix files and
overwrites one stale table with another.

**Two checks that distinguish them, and nothing else does.** The job log's
`flexrv: N v4 files, M fixed, K v7` line -- `K` must equal the ancestry's job
count -- and the table's own `source` column, which must be `v7` for every binary
trait. Note both when a summary is built; after the fact the file is
indistinguishable.

### 14d. A long dsub campaign stalls on accumulated Batch job records (2026-09-30)

dsub never deletes the Batch job it creates, so a production campaign accumulates
records until submissions stop working. **Measured on 2026-09-30: 5,310 Batch jobs
in `us-central1`**, of which 3,524 were pre-today SUCCEEDED, after a day of ~570
reruns plus ~570 splice controls plus summaries.

**The symptom is a HANG, not an error.** `submit()` blocked with no output and the
notebook's interrupt could not clear it, because the block was inside
`subprocess.run` on a `create_job` the client was retrying. `gcloud batch jobs
list` hung in a terminal too -- it pages the whole registry. One resubmission was
silently lost this way and had to be re-issued.

**Diagnosis, from a terminal (the kernel may be wedged):**

```sh
ps aux | grep -E "dsub|gcloud" | grep -v grep | head
gcloud batch jobs list --location us-central1 --format='value(name)' | wc -l
```

**The purge. Deleting a Batch job removes only the job RECORD** -- delivered
outputs in the bucket and the dsub logs under `dsub/logs/` are untouched. Keep the
CURRENT day in every state, because `dstat --full` on those records is the only
history of the day's diagnostics and deleting them makes it unrecoverable.

```sh
nohup bash -c '
  gcloud batch jobs list --location us-central1 \
    --filter="status.state=SUCCEEDED AND createTime<\"YYYY-MM-DDT00:00:00Z\"" \
    --format="value(name)" > /tmp/purge.list
  echo "to delete: $(wc -l < /tmp/purge.list)"
  xargs -n1 -P8 -I{} gcloud batch jobs delete {} --location us-central1 --quiet \
    < /tmp/purge.list
  echo DONE
' > /tmp/purge.log 2>&1 &
```

Detached, because it is one API call per job and takes tens of minutes. Repeat
with `status.state=FAILED` and `CANCELLED` for more headroom.

**Two other things the same evening are probably the same cause**, not
coincidences: repeated `VMs not functioning within the time window 1080 seconds`
followed by requeues, and `CODE_GCE_ZONE_RESOURCE_POOL_EXHAUSTED` on jobs that had
scheduled cleanly earlier in the day. **So purge BEFORE a large campaign, not
after it stalls.**

**REPLICATED on afr (2026-09-30, same evening).** `check7-afr-asth-b0` and
`check7-afr-vte-b0` refuse with the identical signature -- 168 unexplained
`.singleAssoc` rows, `.txt.gz` passing, one gene (`ENSG00000277224`, a DIFFERENT
gene from amr's) and the same `AM_cbrt0.2 x MAC25/MAC50/MAC100` family, whose
pseudo-marker's own count is at or below `--max_MAC_for_ER` (written `<20`) while
its `AM_unw_*` twin's is far larger. So this is **NOT amr-specific**: four
prefixes, two ancestries, two genes, one mechanism, and it will recur on any
ancestry whose genes have rare variants largely unscored by AlphaMissense. Both
afr and amr therefore land at 122 of 124 spliced, with two traits each left
mixed-provenance until the awk change is made. That raises the value of making it:
it is one rule, it unblocks four prefixes now, and it removes a recurring refusal
from every future rollout.

## 15. Timing LEDGER #179 on one chromosome (2026-10-08): old image vs new, same job

For each of a few v4 FlexRV records -- two binary traits of different prevalence
and one quantitative, since only the binary path has saddlepoint and exact-test
cells -- two jobs identical except the image, on one chromosome; then the
durations from dstat and a byte comparison of every output file. The expected
result is "identical" on every file and a shorter wall on the new image; a
different byte is a bug, not a re-baseline. Nothing below prints a count of
people: durations, file sizes and "identical"/"DIFFERENT" only.

```python
# cell 1: submit the pairs. Each v4name is a RUNS key, fxv4-<anc>-<tag(pheno)>-b<i>; the record
# carries anc, pheno, trait, masks, geno and the GRM, so binary and quantitative need no special-casing.
IMG9  = "gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-a08fb295"     # current production
IMG13 = "gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-8bffa0ca"     # LEDGER #179, both sites; built and pushed 2026-10-08, P1 byte-identical in the container
BENCH = ["fxv4-afr-asth-b0", "fxv4-afr-vte-b0", "fxv4-afr-<QUANT-TAG>-b0"]   # two binary, one quantitative
CHRN = "20"
def bench_names(v4name):
    r = RUNS[v4name]; return r, {tag: f"bench179-{r['anc']}-{tag}-{v4name.split('-')[2]}" for tag in ("old", "new")}
for v4name in BENCH:
    r, names = bench_names(v4name); anc = r["anc"]
    env = {k: v for k, v in r["env"].items() if k != "REGIONS"}      # the whole chromosome, no region list
    assert r["geno"] == "pgen", (v4name, r["geno"])
    for tag, img in (("old", IMG9), ("new", IMG13)):
        submit(name=names[tag], image=img, masks=r["masks"], fasttest=r["fasttest"],
               outdir=f"{anc}/bench179_{tag}", anc=anc, pheno=r["pheno"], trait=r["trait"], firth=r["firth"], env=env,
               geno=r["geno"], chrn=CHRN, plk=cod(int(CHRN), anc), grp=am_grp(int(CHRN)), testtype="group",
               pooled=True, sidecars=True, sparse_grm=r["sparse_grm"])
        print(names[tag], r["trait"], RUNS[names[tag]]["geno"])      # geno must read pgen (section 2)
```

```python
# cell 2: durations, once all have finished. Every job includes the same ~60 s of localisation.
import json, subprocess, datetime as dt
import re
AGE = "3h"      # --age: the name filter is applied AFTER the listing; keep it as short as the jobs' age allows
def job_span(name):
    """(status, dstat wall seconds, SAIGE_SECONDS from the job log): the wrapper prints T_SAIGE_START/END
    and SAIGE_SECONDS= around the SAIGE calls, so the third number excludes provisioning, localisation
    and delivery -- the ratio to report. dstat's span is kept as the all-in figure."""
    # GOOGLE_CLOUD_PROJECT (not GOOGLE_PROJECT: empty in a bash -lc shell, and dstat then returns []);
    # --users/--status "*": dstat lists the current user's RUNNING jobs by default; --age: the name
    # filter is applied after the listing, and this project's Batch records make a wide one slow (14d)
    out = subprocess.run(["bash", "-lc", 'dstat --provider google-batch --project "${GOOGLE_CLOUD_PROJECT}" '
                          f"--users '*' --status '*' --age {AGE} --names {name} --full --format json"],
                         capture_output=True, text=True, timeout=120).stdout   # a wide listing can hang the cell (14d);
                         # the kernel cannot interrupt the child -- `pkill -f "dstat --provider"` from a Terminal does
    j = json.loads(out)[0]
    t = lambda k: dt.datetime.fromisoformat(j[k].replace("Z", "+00:00"))
    log = subprocess.run(["gcloud", "storage", "cat", j["logging"]], capture_output=True, text=True).stdout
    m = re.search(r"SAIGE_SECONDS=(\d+)", log)     # no ^ anchor: Batch prefixes every stdout line with [batch_task_logs]...
    return j["status"], (t("end-time") - t("start-time")).total_seconds(), (int(m.group(1)) if m else float("nan"))
print(f"{'record':26s} {'trait':13s} {'old wall':>9s} {'new wall':>9s} {'old SAIGE':>10s} {'new SAIGE':>10s}  ratio(SAIGE)")
for v4name in BENCH:
    r, names = bench_names(v4name)
    (so, ow, os_), (sn, nw, ns) = job_span(names["old"]), job_span(names["new"])
    print(f"{v4name:26s} {r['trait']:13s} {ow:9.0f} {nw:9.0f} {os_:10.0f} {ns:10.0f}  {os_ / ns:5.2f}   [{so}/{sn}]")
```

```python
# cell 3: byte comparison, file by file, nothing displayed but the verdict and sizes. The outputs are
# named after the benchmark JOB (<job name>.<Pheno>.<anc>.chr<N>.<suffix>), not the v4 record. md5_of
# refuses an empty read: the first version swallowed a wrong path and declared two empty strings identical.
import hashlib
def md5_of(path):
    p = subprocess.run(["gcloud", "storage", "cat", path], capture_output=True)
    if p.returncode != 0 or len(p.stdout) == 0:
        raise FileNotFoundError(f"{path}: {p.stderr.decode()[-200:].strip() or 'empty'}")
    return hashlib.md5(p.stdout).hexdigest(), len(p.stdout)
SFX = (".txt.gz", ".singleAssoc.txt.gz", ".pooledTests.txt.gz", ".markerList.txt.gz",
       ".exactByWeight.txt.gz", ".exactByWeightAboveCap.txt.gz", ".spaFallbacks.txt.gz",
       ".spaPinned.txt.gz", ".stretchGate.txt.gz", ".skatFailures.txt.gz", ".skatoMethod.txt.gz", ".chrStatus.txt.gz")
# not .spaDump.tgz: a tar carries mtimes, so its bytes differ even when its contents do not
def bench_path(v4name, tag, sfx):
    r, names = bench_names(v4name)
    return f"{BUCKET}/SAIGE_SLIM_TESTING/{r['anc']}/bench179_{tag}/{names[tag]}.{r['pheno']}.{r['anc']}.chr{CHRN}{sfx}"
for v4name in BENCH:
    r, _ = bench_names(v4name)
    print("==", v4name, r["trait"])
    for sfx in SFX:
        try:
            (h_old, n_old), (h_new, n_new) = md5_of(bench_path(v4name, "old", sfx)), md5_of(bench_path(v4name, "new", sfx))
        except FileNotFoundError as e:
            print(f"  {sfx:32s} MISSING: {e}"); continue
        print(f"  {sfx:32s} {'identical' if h_old == h_new else 'DIFFERENT'}  {n_old} / {n_new} bytes")
```

Gzip embeds no timestamp from SAIGE's writer, but if a `.txt.gz` pair differs
while `zcat | md5sum` agrees, the compressor's header is the difference, not
the rows: check with `gcloud storage cat <path> | zcat | md5sum` on both before
reading a DIFFERENT as a defect.

```python
# cell 4: the same records once more on IMG13 with --is_single_in_groupTest OFF. FlexRV's single-variant
# rows duplicate the SKAT-O run's (real markers) plus one pseudo-marker row per cell that only the splice
# control ever read, so a complete run on one image does not need them. Needs the wrapper patch below
# this section (SINGLE_IN_GROUP, default TRUE: the flag, the OUT_FILE_SINGLE_ASSOC check, and a one-line
# placeholder at delivery because dsub refuses to delocalize a declared output that does not exist).
for v4name in BENCH:
    r, names = bench_names(v4name); anc = r["anc"]
    env = {k: v for k, v in r["env"].items() if k != "REGIONS"}; env["SINGLE_IN_GROUP"] = "FALSE"
    nm = names["new"].replace("-new-", "-nosingle-")
    submit(name=nm, image=IMG13, masks=r["masks"], fasttest=r["fasttest"],
           outdir=f"{anc}/bench179_nosingle", anc=anc, pheno=r["pheno"], trait=r["trait"], firth=r["firth"], env=env,
           geno=r["geno"], chrn=CHRN, plk=cod(int(CHRN), anc), grp=am_grp(int(CHRN)), testtype="group",
           pooled=True, sidecars=True, sparse_grm=r["sparse_grm"])
    print(nm, r["trait"], RUNS[nm]["geno"])
```

```python
# cell 5: the three-way SAIGE seconds and the nosingle-vs-new byte check. Gene-level files must be identical:
# the flag changes nothing a region test reads. .singleAssoc.txt.gz must be the wrapper's one-line
# PLACEHOLDER; the exact-by-weight sidecars may lose their "ultra_rare" rows (that pass runs only with
# the flag on) -- "burden" rows stay.
print(f"{'record':26s} {'trait':13s} {'old':>7s} {'new':>7s} {'nosingle':>9s}   old/new  old/nosingle   (SAIGE seconds)")
for v4name in BENCH:
    r, names = bench_names(v4name)
    (_, _, old), (_, _, new), (_, _, nos) = (job_span(names["old"]), job_span(names["new"]),
                                             job_span(names["new"].replace("-new-", "-nosingle-")))
    print(f"{v4name:26s} {r['trait']:13s} {old:7.0f} {new:7.0f} {nos:9.0f}   {old / new:5.2f}     {old / nos:5.2f}")
for v4name in BENCH:
    r, names = bench_names(v4name); anc = r["anc"]
    nos = names["new"].replace("-new-", "-nosingle-")
    nos_path = lambda sfx: f"{BUCKET}/SAIGE_SLIM_TESTING/{anc}/bench179_nosingle/{nos}.{r['pheno']}.{anc}.chr{CHRN}{sfx}"
    print("==", v4name, r["trait"])
    for sfx in (".txt.gz", ".pooledTests.txt.gz", ".markerList.txt.gz"):
        (ha, na), (hb, nb) = md5_of(bench_path(v4name, "new", sfx)), md5_of(nos_path(sfx))
        print(f"  {sfx:24s} {'identical' if ha == hb else 'DIFFERENT'}  {na} / {nb} bytes")
    single = nos_path(".singleAssoc.txt.gz")
    head = subprocess.run(["bash", "-lc", f"gcloud storage cat {single} | zcat | head -c 200"], capture_output=True, text=True).stdout
    ok = head.startswith("# no singleAssoc output")
    print(f"  .singleAssoc.txt.gz      {'placeholder, as it must be' if ok else 'REAL ROWS -- the flag did not reach SAIGE'}")
```

Expected from sim50k: old/new ~1.5, old/nosingle ~2.8 on SAIGE's own seconds
(the wall includes provisioning, localisation and delivery on every side).

The wrapper patch for cell 4 (`saige_step2_brava_main_original.sh`; every
existing job is unchanged because the default is TRUE and the TRUE branches are
the old code):

```bash
# inputs block
SINGLE_IN_GROUP="${SINGLE_IN_GROUP:-TRUE}"            # TRUE | FALSE: --is_single_in_groupTest (group tests)

# group-test checks: replace the unconditional OUT_FILE_SINGLE_ASSOC requirement
  case "${SINGLE_IN_GROUP}" in TRUE|FALSE) ;; *) die "SINGLE_IN_GROUP must be TRUE or FALSE (got '${SINGLE_IN_GROUP}')" ;; esac
  if [[ -z ${MULTI_MODE} && ${SINGLE_IN_GROUP} == "TRUE" ]]; then
    [[ -n ${OUT_FILE_SINGLE_ASSOC} ]] || die "OUT_FILE_SINGLE_ASSOC not set for a group test"
  fi

# the command's group block
    --is_single_in_groupTest="${SINGLE_IN_GROUP}"

# delivery: replace the two singleAssoc lines
  if [[ ${SINGLE_IN_GROUP} == "TRUE" ]]; then
    catparts ".singleAssoc.txt" || die "no chromosome produced a singleAssoc file"
    deliver "${OUTBASE}.singleAssoc.txt" "${OUT_FILE_SINGLE_ASSOC}"
  elif [[ -n ${OUT_FILE_SINGLE_ASSOC} ]]; then
    # declared output, so it must exist for dsub to delocalize (same reason as .pooledTests)
    printf '# no singleAssoc output: SINGLE_IN_GROUP=FALSE (variant results come from the SKAT-O run)\n' \
      > "${OUTBASE}.singleAssoc.txt"
    deliver "${OUTBASE}.singleAssoc.txt" "${OUT_FILE_SINGLE_ASSOC}"
  fi
```

A job that stays RUNNING with an EMPTY log is not diagnosed by the log: dsub
copies the task's output from Cloud Logging into the GCS log when the task
ends, so an empty file means only "not finished". Read the record's `events`
instead (`dstat ... --full --format json`, keys `events`, `last-update`,
`status-message`): a `scheduled` event followed by nothing, or by a
`CODE_GCE_ZONE_RESOURCE_POOL_EXHAUSTED` message, is the zone queue (seen
2026-10-08 on one of three bench179 jobs, us-central1-c, three failed VM
creations in 2.5 min; Batch keeps retrying). Cancel and resubmit if it matters;
a fresh submission can land elsewhere.

## 16. The spliceAI rule test (2026-10-12): old and new rules in ONE SKAT-O run

`~/Repositories/brava-reannotation/out/aou/spliceai_split_test/aou.chr<c>.spliceai_split_test.saige_group.txt.gz`
(built by `build/split_spliceai_test.py`): the production variant set minus the gnomAD popmax > 0.01
drops, same `chr:pos:ref:alt` IDs, pan-ancestry, NO score line (SKAT-O only), with the old damaging
set split exactly three ways -- `damaging_missense` (REVEL/CADD/LOFTEE LC), `spliceai_gt05`
(DS >= 0.5), `spliceai_02_05` (DS 0.2-0.5). `other_missense` and `non_coding` here LACK the DS
0.2-0.5 variants (as under the old rule). On chr20: 37,429 / 6,770 / 15,878; against production's
own labels the old set differs by ~5-9k pairs (dbNSFP 5.4a, the multi-gene spliceAI fix), so "old"
means the old RULE inside the new pipeline. `docs/spliceai.md` in that repo has the calibration
behind 0.5.

Masks (owner, 2026-10-12), SKAT-O at the production MAF bins, ten:
`pLoF,synonymous,damaging_missense,damaging_missense:spliceai_gt05,damaging_missense:spliceai_gt05:spliceai_02_05,pLoF:damaging_missense:spliceai_gt05,spliceai_gt05,spliceai_02_05,other_missense,pLoF:damaging_missense:spliceai_gt05:spliceai_02_05:other_missense:synonymous`.
What each is for: `pLoF` and the all-in mask (the last) are unchanged by the rule and must
reproduce production up to the pipeline differences (the positive control on the file and the run);
`synonymous` is the null beside the power claim; `damaging_missense:spliceai_gt05` (new rule) against
`damaging_missense:spliceai_gt05:spliceai_02_05` (old rule) is the claim itself, gene by gene with
only the moved variants differing; `spliceai_02_05` alone tests the calibration -- near null if DS
0.2-0.5 is only Supporting-level evidence, hits in known genes if the cutoff discards recall;
`other_missense` is the old-rule version (without the DS 0.2-0.5 variants).

```sh
# upload (laptop -> bucket; uploads are free). 22 autosomes, ~380 MB.
BUCKET=gs://<the workspace bucket>
gcloud storage cp ~/Repositories/brava-reannotation/out/aou/spliceai_split_test/aou.chr{1..22}.spliceai_split_test.saige_group.txt.gz \
  "$BUCKET/SAIGE_SLIM_TESTING/annotations_reannot_v1/"
gcloud storage ls "$BUCKET/SAIGE_SLIM_TESTING/annotations_reannot_v1/" | wc -l     # 22
```

```python
# submit: copy the production SKAT-O record (skv4-<anc>-<tag>-b<i>) and change ONLY the group file,
# the masks and the output directory. Pilot = afr chr20, the three benchmark traits; then all.
IMG13 = "gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-8bffa0ca"
SPLIT_MASKS = ("pLoF,synonymous,damaging_missense,damaging_missense:spliceai_gt05,"
               "damaging_missense:spliceai_gt05:spliceai_02_05,pLoF:damaging_missense:spliceai_gt05,"
               "spliceai_gt05,spliceai_02_05,other_missense,"
               "pLoF:damaging_missense:spliceai_gt05:spliceai_02_05:other_missense:synonymous")
split_grp = lambda c: f"{BUCKET}/SAIGE_SLIM_TESTING/annotations_reannot_v1/aou.chr{c}.spliceai_split_test.saige_group.txt.gz"
PILOT = ["skv4-afr-asth-b0", "skv4-afr-vte-b0", "skv4-afr-bmi-b0"]; CHRN = "20"
for v4name in PILOT:
    r = RUNS[v4name]; anc = r["anc"]
    assert r["geno"] == "pgen" and "RCORR" not in r["env"], v4name      # a SKAT-O record: r.corr at the build default
    env = {k: v for k, v in r["env"].items() if k != "REGIONS"}
    nm = v4name.replace("skv4", "sai05")
    submit(name=nm, image=IMG13, masks=SPLIT_MASKS, fasttest=r["fasttest"],
           outdir=f"{anc}/skato_sai05", anc=anc, pheno=r["pheno"], trait=r["trait"], firth=r["firth"], env=env,
           geno=r["geno"], chrn=CHRN, plk=cod(int(CHRN), anc), grp=split_grp(int(CHRN)), testtype="group",
           pooled=True, sidecars=True, sparse_grm=r["sparse_grm"])
    print(nm, r["trait"], RUNS[nm]["geno"])
```

Reading the pilot: `pLoF` and the all-in mask against the production SKAT-O run for the same gene
(agreement up to the pipeline differences bounds what the rule itself moves); `synonymous` on the
null; then `damaging_missense:spliceai_gt05` against `damaging_missense:spliceai_gt05:spliceai_02_05`
gene by gene, and `spliceai_02_05` alone. Afterwards, per mask: genes past exome-wide significance,
replication of the known gene-trait pairs, and the inflation beside each (section 7c's calibration
column takes these rows), so a gain is not read off an inflated tail.
Full run: the same cell over every `skv4-*` record and `chrs=range(1, 23)` with `grp=split_grp`
(the function, as section 2 does for multi-chromosome jobs). Single-in-group stays as the record had
it: this arm's single-variant file IS the variant results for everything downstream.

**Pilot outcome (2026-10-12, afr chr20, asth/vte/bmi):** `pLoF` cells identical to production in every
gene; all-in and `synonymous` 96-99 % identical, the rest shifted <= 0.9 log10 (popmax drops, dbNSFP
relabels); synonymous null 0.048-0.065 at 0.05; new vs old damaging a coin toss gene by gene, nothing
past 2.5e-6 under either, the dropped set's smallest p 1.8e-4 over ~1,600 cells. Plumbing passes;
chr20 x three traits cannot show the rule. The evidence run is below.

```python
# the evidence run: eur, genome-wide, the ten masks, eight quantitative traits with many established
# rare-variant genes. One job per production bin (skv4-eur-<tag>-b0..b3, each its own chromosome
# set, machine type and spot setting), changing ONLY the group file, the masks, the image and the
# output directory. eur has ten bins per trait: 80 jobs, spot. Reading: spliceai_02_05 pooled over all genes x traits on a QQ (null if
# the cutoff is right); the new vs old damaging mask at the Genebass burden genes for these traits;
# exome-wide hits per mask with the inflation beside each; synonymous as the null throughout.
TRAITS = ["ldlc", "hdlc", "tg", "tchol", "alt", "hba1c", "height", "bmi"]
def bins_of(anc, tag):      # eur was binned into MORE than four jobs; take whatever the registry holds, and check it tiles 1-22
    ks = sorted(k for k in RUNS if k.startswith(f"skv4-{anc}-{tag}-b"))
    chrs = [c for k in ks for c in RUNS[k]["chrs"]]
    assert sorted(chrs) == list(range(1, 23)), (tag, ks, sorted(chrs))
    return ks
for tag in TRAITS:
    for v4name in bins_of("eur", tag):
        r = RUNS[v4name]; anc = r["anc"]
        assert r["geno"] == "pgen" and "RCORR" not in r["env"], v4name
        env = {k: v for k, v in r["env"].items() if k != "REGIONS"}
        nm = v4name.replace("skv4", "sai05")
        submit(name=nm, image=IMG13, masks=SPLIT_MASKS, fasttest=r["fasttest"],
               outdir=f"{anc}/skato_sai05", anc=anc, pheno=r["pheno"], trait=r["trait"], firth=r["firth"], env=env,
               geno=r["geno"], chrs=r["chrs"], plk=lambda c, a=anc: cod(c, a), grp=split_grp, testtype="group",
               machine_type=r["machine_type"], spot=r["spot"], pooled=True, sidecars=True, sparse_grm=r["sparse_grm"])
        # 2026-10-12: every b0 (chr1 alone) failed in ~60 s -- the upload had dropped aou.chr1's file, the
        # job log's `gcloud.storage.cp ... matched no objects` named it. Count the listing (22) before submitting,
        # and check SIZES: the workbench's web uploader stalled on the re-upload and left a 9,437,184-byte object
        # (exactly 9 MiB, a chunk boundary) for a 19,128,189-byte file -- a truncated gzip that lists like a
        # file. `gcloud storage cp` from the laptop is the route that works.
        # spot as production ran it (owner, 2026-10-12): preempted jobs are resubmitted from the check below
        print(nm, r["chrs"], RUNS[nm]["geno"])
```

```python
# after the batch settles: which jobs to resubmit. A preempted job and a failed one look alike from
# outside; and a job can report SUCCESS with chromosomes missing (the wrapper exits 0 so the ones that
# worked are delivered) -- only chrStatus says so. Resubmit = the cell above filtered to `redo`.
def status_of(name):
    out = subprocess.run(["bash", "-lc", 'dstat --provider google-batch --project "${GOOGLE_CLOUD_PROJECT}" '
                          f"--users '*' --status '*' --age 2d --names {name} --format json"], capture_output=True, text=True, timeout=120).stdout
    j = json.loads(out); return j[0]["status"] if j else "NOT FOUND"
redo = []
for tag in TRAITS:
    for v4name in bins_of("eur", tag):
        nm = v4name.replace("skv4", "sai05"); st = status_of(nm)
        if st != "SUCCESS": redo.append(nm); print(nm, st); continue
        chs = rows(RUNS[nm]["out_prefix"] + ".chrStatus.txt.gz")
        failed = [l for l in chs if l.startswith("failed")][0].split("\t")[1]
        if failed != "none": redo.append(nm); print(nm, "SUCCESS but chromosomes failed:", failed)
print(len(redo), "to resubmit:", redo)
```

```python
# reading the evidence run (after every bin is in). Pools the ten bins per trait. Prints counts, gene
# IDs and p-values only. Reference genes = exome-wide significant under pLoF in this run: pLoF is
# unchanged by the rule, so the set is independent of the comparison and it is where a missense/splice
# mask has signal to gain. (No local Genebass gene list; results.mt is in the cloud.)
import math
NEW_DM, OLD_DM, DROPPED = "damaging_missense;spliceai_gt05", "damaging_missense;spliceai_gt05;spliceai_02_05", "spliceai_02_05"
SIG = 2.5e-6
def load_trait(tag):
    t = {}
    for v4name in bins_of("eur", tag):
        nm = v4name.replace("skv4", "sai05"); t.update(table(rows(RUNS[nm]["out_prefix"] + ".txt.gz")))
    return t
def minp_by_gene(t, mask):
    d = {}
    for (g, grp, m), p in t.items():
        if grp == mask and p != "NA": d[g] = min(d.get(g, 1.0), float(p))
    return d
def qq_quantiles(ps):   # -log10 observed at the expected 1e-1 .. 1e-5 points, plus the fraction below 0.05
    ps = sorted(ps); n = len(ps); out = []
    for e in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5):
        k = int(e * n); out.append((e, ps[k] if k < n else float("nan")))
    return n, out, sum(p < 0.05 for p in ps) / n
pooled_dropped, pooled_syn = [], []
print(f"{'trait':7s} {'genes':>6s} {'pLoF ref':>8s} | {'new<old':>7s} {'new>old':>7s} | {'hits new':>8s} {'hits old':>8s} {'hits dropped':>12s} {'hits allin':>10s} | ref genes: new p < old p")
for tag in TRAITS:
    t = load_trait(tag)
    pn, po, pd = minp_by_gene(t, NEW_DM), minp_by_gene(t, OLD_DM), minp_by_gene(t, DROPPED)
    plof = minp_by_gene(t, "pLoF"); allin = minp_by_gene(t, [k[1] for k in t if k[1].startswith("pLoF;damaging_missense;spliceai_gt05;spliceai_02_05")][0])
    ref = [g for g, p in plof.items() if p < SIG]
    common = set(pn) & set(po)
    better = sum(pn[g] < po[g] for g in common); worse = sum(pn[g] > po[g] for g in common)
    ref_better = sum(pn.get(g, 1) < po.get(g, 1) for g in ref); ref_n = sum(g in pn and g in po for g in ref)
    print(f"{tag:7s} {len(common):6d} {len(ref):8d} | {better:7d} {worse:7d} | {sum(p < SIG for p in pn.values()):8d} {sum(p < SIG for p in po.values()):8d} "
          f"{sum(p < SIG for p in pd.values()):12d} {sum(p < SIG for p in allin.values()):10d} | {ref_better} of {ref_n}")
    pooled_dropped += [float(p) for (g, grp, m), p in t.items() if grp == DROPPED and p != "NA"]
    pooled_syn += [float(p) for (g, grp, m), p in t.items() if grp == "synonymous" and p != "NA"]
    if ref:
        print("   ref genes (pLoF-significant), new vs old damaging p:",
              [(g, f"{pn.get(g, float('nan')):.1e}", f"{po.get(g, float('nan')):.1e}") for g in sorted(ref, key=lambda g: plof[g])[:8]])
for label, ps in (("spliceai_02_05 pooled", pooled_dropped), ("synonymous pooled", pooled_syn)):
    n, q, f05 = qq_quantiles(ps)
    print(f"{label}: n={n}, frac<0.05={f05:.3f}; observed p at expected " + ", ".join(f"{e:.0e}: {o:.1e}" for e, o in q))
```

Reading: `synonymous pooled` must track its expected points (the null); `spliceai_02_05 pooled`
tracking them says the dropped variants carry no burden signal and the cutoff loses nothing; lifting
off at 1e-3..1e-5 says it discards recall. The hits columns and the ref-gene comparison are the power
claim; the all-in column is the tie to production. A gain in hits with an inflated synonymous line is
not a gain.

**Evidence-run outcome (2026-10-13, eur genome-wide, 8 quantitative traits, 80 bins + 2 resubmitted
on demand after spot preemption; ~18k genes a trait):**

| trait | pLoF ref genes | new<old / new>old (all genes) | hits new | hits old | hits dropped | hits all-in | ref genes new<old |
|---|---|---|---|---|---|---|---|
| ldlc | 4 | 8022 / 8248 | 2 | 2 | 0 | 5 | 4 of 4 |
| hdlc | 10 | 8324 / 8102 | 8 | 6 | 0 | 13 | 4 of 10 |
| tg | 7 | 8321 / 8110 | 6 | 5 | 0 | 12 | 3 of 7 |
| tchol | 6 | 8134 / 8307 | 4 | 5 | 1 | 8 | 6 of 6 |
| alt | 2 | 8294 / 8211 | 1 | 1 | 0 | 2 | 1 of 2 |
| hba1c | 1 | 7897 / 8270 | 2 | 1 | 0 | 1 | 0 of 1 |
| height | 25 | 8565 / 8225 | 25 | 24 | 1 | 49 | 14 of 25 |
| bmi | 7 | 8396 / 8323 | 2 | 1 | 0 | 2 | 2 of 7 |
| t2diab (binary, 2026-10-13) | 4 | 8235 / 8462 | 3 | 2 | 0 | 1 | 3 of 4 |

Pooled QQ (observed p at the expected 1e-1 / 1e-2 / 1e-3 / 1e-4 / 1e-5 points; frac < 0.05):
`spliceai_02_05` n=402,419: 9.7e-2 / 8.5e-3 / 8.1e-4 / 5.7e-5 / 5.9e-6; 0.052 -- ON THE NULL.
`synonymous` n=441,731: 9.5e-2 / 7.7e-3 / 5.3e-4 / 2.9e-5 / 1.7e-7; 0.055 (the usual large-gene
synonymous tail; heavier than the dropped set's). New damaging n=435,979: 9.0e-2 / 6.3e-3 / 2.3e-4 /
1.2e-11 / 9.4e-41; old damaging n=439,501: 9.0e-2 / 6.2e-3 / 2.6e-4 / 2.1e-11 / 7.1e-41 -- the same tail.

t2diab, the binary check (one trait, ten bins on demand): `spliceai_02_05` n=51,098, 0.046 below 0.05,
QQ 1.1e-1 / 9.3e-3 / 8.4e-4 / 4.2e-5 / 3.5e-5 (at or above expected); `synonymous` 0.054; per MAF bin
both masks 0.044-0.055 at 0.05 and 0.0011-0.0017 at 0.001, so the binary tail's few-case inflation is
the same under both rules; new/old damaging tails 8.8e-13 / 1.0e-12 at the 1e-5 point; HNF1A
2.2e-4 -> 8.7e-5.

Reading: the DS 0.2-0.5 variants carry no burden signal in aggregate (2 hits in 8 traits against 50);
the new mask finds 50 exome-wide hits against the old rule's 45, new >= old on 7 of 8 traits, with
identical tail behaviour and the same synonymous null; at the 62 pLoF-reference genes the new mask is
better at 34 -- gains where the rule cleans (APOC3 3.1e-4 -> 4.2e-7 for HDLC, LDLR, LCAT, ACAN),
losses where some 0.2-0.5 splice variants were real (CD36 2.2e-8 -> 9.0e-7, CETP, one height gene
3.6e-7 -> 3.3e-5). **DECIDED (owner, 2026-10-13, on eight quantitative traits and one binary): the spliceAI cutoff moves to 0.5**; production runs
the seven-mask set on the release `new_rules` file:
`pLoF,damaging_missense,spliceai_gt05,pLoF:damaging_missense:spliceai_gt05,synonymous,other_missense,pLoF:damaging_missense:spliceai_gt05:other_missense:synonymous`
(the all-in mask under the new rule carries the DS 0.2-0.5 variants inside `other_missense`/`non_coding`,
so on the release file it is spelt without `spliceai_02_05`). eur bin times on IMG13, ten masks, chromosomes 3+8 (`-b5`), SAIGE seconds from the job log:
quantitative (ldlc) 972 s, binary (t2diab) 3,416 s -- 3.5x, the per-cell saddlepoint of the binary
SKAT-O path. So a production rerun of ~60 binary traits is ~600 bin-hours on n2-highmem-2 (~$60
on demand), the ~12 quantitative traits a sixth of that; spot lost 2 of 80 quantitative bins tonight.

## 17. SpliceAI v1 (DS >= 0.5) against SpliceAI2 (>= 0.25) in ONE SKAT-O run (2026-10-14)

`~/Repositories/brava-reannotation/out_spliceai2/aou/sai_versions_test/aou.chr<c>.sai_versions_test.saige_group.txt.gz`
(`build/sai_versions_test.py`): every entry carries its SpliceAI2-build label, except the splice category,
which is split by which version calls it -- `sai_both`, `sai_v1only` (v1 >= 0.5 only), `sai_v2only`
(SpliceAI2 >= 0.25 only). The builder refuses to write if `pLoF` or `damaging_missense` differ between the
builds (neither involves SpliceAI). Genome-wide: v1 27,326, v2 36,990, both 18,884 (Jaccard 0.42); chr20
2,045 / 4,724 / 4,581. chr1 ships as two parts (`.part_aa` 15,728,640 + `.part_ab` 3,344,665 bytes) because
the workbench uploader cannot take the 19 MB file; `cat` of the parts has the full file's md5.

Masks, twelve: the three unchanged controls, the two rules, each version's splice set, the three pieces
of the split (the disagreement sets `sai_v1only` / `sai_v2only` are the decisive ones), `other_missense`,
and the all-in mask (under both rules at once, since it holds every splice piece):
`pLoF,synonymous,damaging_missense,other_missense,damaging_missense:sai_both:sai_v1only,damaging_missense:sai_both:sai_v2only,sai_both:sai_v1only,sai_both:sai_v2only,sai_both,sai_v1only,sai_v2only,pLoF:damaging_missense:sai_both:sai_v1only:sai_v2only:other_missense:synonymous`

Upload (workbench web uploader, into `SAIGE_SLIM_TESTING/annotations_reannot_v2/`): chr2..22 as they are,
chr1 as its two parts; then in the notebook compose chr1 and md5 everything against the laptop list:

```python
import hashlib, subprocess
D = f"{BUCKET}/SAIGE_SLIM_TESTING/annotations_reannot_v2"
# chr1 from its parts, server side; each part's size checked first (compose would join a truncated part)
parts = [f"{D}/aou.chr1.sai_versions_test.saige_group.txt.gz.part_{s}" for s in ("aa", "ab")]
for p, n in zip(parts, (15728640, 3344665)):
    got = len(subprocess.run(["gcloud", "storage", "cat", p], capture_output=True).stdout); assert got == n, (p, got, n)
dest = f"{D}/aou.chr1.sai_versions_test.saige_group.txt.gz"
print(subprocess.run(["gcloud", "storage", "objects", "compose", *parts, dest], capture_output=True, text=True).stderr[-200:])
REF = """ea5c3e2c33eb156220cee2110e80a291 19073305 aou.chr1.sai_versions_test.saige_group.txt.gz
b23baa838835bd7392cad3bcb130c228 13025514 aou.chr2.sai_versions_test.saige_group.txt.gz
4f2a88ce097656bf57f90b1e75d698d0 11184831 aou.chr3.sai_versions_test.saige_group.txt.gz
088b1acc793aa61e89304265ceb2040d 7191211 aou.chr4.sai_versions_test.saige_group.txt.gz
eed8eb71d3862a618ae3e3631514222c 11298945 aou.chr5.sai_versions_test.saige_group.txt.gz
3cec5daa349027cf8f2c1a6c58c7e2a6 9628159 aou.chr6.sai_versions_test.saige_group.txt.gz
50f671281b1c71c383d1e558bca693a5 8855502 aou.chr7.sai_versions_test.saige_group.txt.gz
e5fb486253ad5ba3163735fd9c812df5 6717005 aou.chr8.sai_versions_test.saige_group.txt.gz
bcb981f3bd67e7c25a6742ddb03ecce3 7763187 aou.chr9.sai_versions_test.saige_group.txt.gz
eae14fc1419da01294fe5878b4e85f1b 7201901 aou.chr10.sai_versions_test.saige_group.txt.gz
93de48441a84b320beddedecf9626ee7 11754107 aou.chr11.sai_versions_test.saige_group.txt.gz
6f4d2c3a5c044a477203df2bb59ce160 10345765 aou.chr12.sai_versions_test.saige_group.txt.gz
547357e042edcf091d40b990f8f3811f 3437785 aou.chr13.sai_versions_test.saige_group.txt.gz
021fdcb9813adb8de629e76b7f9505ac 6492259 aou.chr14.sai_versions_test.saige_group.txt.gz
069287d791149c00aecb0eea5d528185 6910637 aou.chr15.sai_versions_test.saige_group.txt.gz
9a42c5f1267cb038333c7958816299fb 9573142 aou.chr16.sai_versions_test.saige_group.txt.gz
74d2ebaafa65c978f9387707af237a0e 11561726 aou.chr17.sai_versions_test.saige_group.txt.gz
bfc6ccf83f973b0c939a4a0f915b5bc7 3032049 aou.chr18.sai_versions_test.saige_group.txt.gz
79f064c0ee55c1d1f61178976bda44d8 13363704 aou.chr19.sai_versions_test.saige_group.txt.gz
c1b07d370f1c0223e0a50e21b707de24 4853842 aou.chr20.sai_versions_test.saige_group.txt.gz
812dd545da0df3eb43300209b963864d 1838382 aou.chr21.sai_versions_test.saige_group.txt.gz
b1c0cc3c11413d42ac6db047464a868f 4152616 aou.chr22.sai_versions_test.saige_group.txt.gz"""
bad = []
for l in REF.splitlines():
    h, n, f = l.split(); data = subprocess.run(["gcloud", "storage", "cat", f"{D}/{f}"], capture_output=True).stdout
    ok = hashlib.md5(data).hexdigest() == h and len(data) == int(n)
    print(f"{f:50s} {len(data):>10d} / {n:<10s} {'ok' if ok else 'MISMATCH'}", flush=True); bad += [] if ok else [f]
print(len(bad), "bad:", bad)
if not bad:
    for p in parts: subprocess.run(["gcloud", "storage", "rm", p])
    print("chr1 parts removed")
```

```python
# submit: eur, the eight quantitative traits + T2D, every bin, same image; only the group file, masks and
# output directory differ from section 16. (Spot: 2 of 80 lost last time; the check cell lists them.)
SAI2_MASKS = ("pLoF,synonymous,damaging_missense,other_missense,"
              "damaging_missense:sai_both:sai_v1only,damaging_missense:sai_both:sai_v2only,"
              "sai_both:sai_v1only,sai_both:sai_v2only,sai_both,sai_v1only,sai_v2only,"
              "pLoF:damaging_missense:sai_both:sai_v1only:sai_v2only:other_missense:synonymous")
sai2_grp = lambda c: f"{BUCKET}/SAIGE_SLIM_TESTING/annotations_reannot_v2/aou.chr{c}.sai_versions_test.saige_group.txt.gz"
TRAITS2 = ["ldlc", "hdlc", "tg", "tchol", "alt", "hba1c", "height", "bmi", "t2diab"]
for tag in TRAITS2:
    for v4name in bins_of("eur", tag):
        r = RUNS[v4name]; anc = r["anc"]
        assert r["geno"] == "pgen" and "RCORR" not in r["env"], v4name
        env = {k: v for k, v in r["env"].items() if k != "REGIONS"}
        nm = v4name.replace("skv4", "sai2v")
        submit(name=nm, image=IMG13, masks=SAI2_MASKS, fasttest=r["fasttest"],
               outdir=f"{anc}/skato_sai2v", anc=anc, pheno=r["pheno"], trait=r["trait"], firth=r["firth"], env=env,
               geno=r["geno"], chrs=r["chrs"], plk=lambda c, a=anc: cod(c, a), grp=sai2_grp, testtype="group",
               machine_type=r["machine_type"], spot=r["spot"], pooled=True, sidecars=True, sparse_grm=r["sparse_grm"])
        print(nm, r["chrs"], RUNS[nm]["geno"])
```

The delivery check is section 16's listing cell with `skato_sai2v` and `sai2v`; the reading cell is below.

```python
# reading: v1 against v2, per trait and pooled. Decisive lines: sai_v1only and sai_v2only pooled -- the
# variants only one version calls; a version whose private set carries signal is finding real splice
# damage the other misses, one whose private set is null is adding noise.
import gzip, math, subprocess
def rows(path):
    data = subprocess.run(["gcloud", "storage", "cat", path], capture_output=True, check=True).stdout
    return gzip.decompress(data).decode().splitlines()
def table(lines):
    h = lines[0].split("\t"); ix = {k: i for i, k in enumerate(h)}; out = {}
    for l in lines[1:]:
        f = l.split("\t")
        if f[ix["Group"]] == "Cauchy": continue
        out[(f[ix["Region"]], f[ix["Group"]], f[ix["max_MAF"]])] = f[ix["Pvalue"]]
    return out
listing = subprocess.run(["gcloud", "storage", "ls", f"{BUCKET}/SAIGE_SLIM_TESTING/eur/skato_sai2v/"], capture_output=True, text=True).stdout.split()
prefix_of = {p.rsplit("/", 1)[-1].split(".")[0]: p[:-len(".chrStatus.txt.gz")] for p in listing if p.endswith(".chrStatus.txt.gz")}
V1, V2 = "damaging_missense;sai_both;sai_v1only", "damaging_missense;sai_both;sai_v2only"
S1, S2, BOTH, ONLY1, ONLY2 = "sai_both;sai_v1only", "sai_both;sai_v2only", "sai_both", "sai_v1only", "sai_v2only"
SIG = 2.5e-6
def minp_by_gene(t, mask):
    d = {}
    for (g, grp, m), p in t.items():
        if grp == mask and p != "NA": d[g] = min(d.get(g, 1.0), float(p))
    return d
def qq(ps):
    ps = sorted(ps); n = len(ps)
    return n, [(e, ps[int(e * n)]) for e in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5) if int(e * n) < n], sum(p < 0.05 for p in ps) / n
pooled = {m: [] for m in (ONLY1, ONLY2, BOTH, S1, S2, V1, V2, "synonymous")}
print(f"{'trait':7s} {'ref':>4s} | {'v1<v2':>6s} {'v1>v2':>6s} | {'hits v1':>7s} {'hits v2':>7s} {'splice v1':>9s} {'splice v2':>9s} {'v1only':>6s} {'v2only':>6s} {'both':>5s} | ref genes v2<v1")
for tag in TRAITS2:
    t = {}
    for v4name in bins_of("eur", tag):
        nm = v4name.replace("skv4", "sai2v"); assert nm in prefix_of, nm; t.update(table(rows(prefix_of[nm] + ".txt.gz")))
    p1, p2, s1, s2, pb, o1, o2, plof = (minp_by_gene(t, m) for m in (V1, V2, S1, S2, BOTH, ONLY1, ONLY2, "pLoF"))
    ref = [g for g, p in plof.items() if p < SIG]; common = set(p1) & set(p2)
    hits = lambda d: sum(p < SIG for p in d.values())
    print(f"{tag:7s} {len(ref):4d} | {sum(p1[g] < p2[g] for g in common):6d} {sum(p1[g] > p2[g] for g in common):6d} | {hits(p1):7d} {hits(p2):7d} {hits(s1):9d} {hits(s2):9d} {hits(o1):6d} {hits(o2):6d} {hits(pb):5d} | "
          f"{sum(g in p1 and g in p2 and p2[g] < p1[g] for g in ref)} of {sum(g in p1 and g in p2 for g in ref)}", flush=True)
    if ref: print("   ref genes, v1 vs v2 damaging p:", [(g, f"{p1.get(g, float('nan')):.1e}", f"{p2.get(g, float('nan')):.1e}") for g in sorted(ref, key=lambda g: plof[g])[:8]])
    for m in pooled: pooled[m] += [float(p) for (g, grp, mm), p in t.items() if grp == m and p != "NA"]
for m, ps in pooled.items():
    n, q, f05 = qq(ps)
    print(f"{m:36s} n={n:7d} frac<0.05={f05:.3f}  observed at expected " + ", ".join(f"{e:.0e}: {o:.1e}" for e, o in q))
```

### 17b. SpliceAI2's cut-off moved to 0.40 and the run was redone (2026-10-15)

The reannotation side recalibrated SpliceAI2 against ACMG (`brava-reannotation`
`spliceai.md`, commits `717f86b`/`8e32fe7`, 2026-10-14): 0.25 is the README's
"Supporting" line, 0.40 is "Moderate", and the rule is now `spliceai2_gt040`. The
split-test files were rebuilt at 0.40 (same three labels, `sai_both` /
`sai_v1only` / `sai_v2only`; chr21: 1,466 / 1,311 / 692 variants, Jaccard 0.42
between the versions' calls), uploaded to `SAIGE_SLIM_TESTING/annotations_reannot_v3/`
(chr1 composed from two parts, every object md5-checked against
`out_spliceai2/aou/sai_versions_test.md5.txt`), and the section-17 cells were
rerun as job names `sai2v040-...`, output `eur/skato_sai2v040`, on IMG14
(`LOG10P` on, so every file carries the `log10_Pvalue*` columns; the reading
cell reads by name). The generated cells are
`brava-reannotation/workbench/sai2v040_cells.md`. **The 0.25 run (section 17,
`skato_sai2v`) is superseded**: its six read traits were a near draw with v2
slightly ahead, and its height/bmi/t2diab outputs are not read. The decision --
v1 at 0.5 against v2 at 0.40 -- is taken on this run's reading cell over all
nine traits, the private sets (`sai_v1only`, `sai_v2only`) pooled being the
decisive lines.

**Outcome (2026-10-15, all nine traits read; DECIDED: SpliceAI2 at 0.40).** The decisive
lines are null for BOTH versions: over the eight non-height traits `sai_v1only`
(294,717 cells) and `sai_v2only` (198,010) track the synonymous line at 1e-2..1e-4
(v1's 2.7e-7 at the 1e-5 quantile is its single hdlc hit, 0.7 such cells expected by
chance at that n; v2's private set has no hit), and on height, whose synonymous control
is itself far off the null (min p 1.1e-15), `sai_v1only` tracks synonymous with 2 hits
while `sai_v2only` sits BELOW it with 1. Neither version's exclusive calls carry signal
these traits can see, so the choice is a tie-break, and it is close either way:
damaging + splice hits 51 (v1) vs 54 (v2) over nine traits (height 25 vs 27); gene for
gene, v2's p is the smaller in a majority for all nine traits (by 0.5-3 %); the pooled
non-height tail is sharper under v2 at every quantile (1e-4: 4.2e-7 vs 5.9e-7; 1e-5:
1.2e-35 vs 1.5e-34), height's is mixed (1e-4 v1, 1e-5 v2); the pLoF-reference genes
favour v1, 38 of 65 (APOB, CETP, GIGYF1 sharpen under v1; APOC3, CD36, LCAT, ABCB4
under v2). Verdict per the rule stated before the numbers: a draw on the private sets
goes to the hit count and the calibrated cut-off, i.e. v2. The cost of being wrong is
about one hit either way. One observation to keep: even `sai_both` is barely above
synonymous outside height (4.2e-5 vs 3.6e-5 at 1e-4) and the splice-only masks have
zero non-height hits -- splice annotation earns its place inside the combined damaging
mask, not on its own.

**Production mask set under the decision** (release files
`brava-reannotation/out_spliceai2/aou/aou.chr<c>.new_rules.saige_group.txt.gz`, labels
`pLoF damaging_missense spliceai2_gt040 other_missense synonymous non_coding`; NOT yet
uploaded):
`pLoF,damaging_missense,pLoF:damaging_missense,pLoF:damaging_missense:spliceai2_gt040,synonymous,other_missense,pLoF:damaging_missense:spliceai2_gt040:other_missense:synonymous`
(owner, 2026-10-15: the with/without pair is `pLoF:damaging_missense` against
`pLoF:damaging_missense:spliceai2_gt040`; a splice-only mask was null on its own in every run
and is not spent).
`am_moderate` is its own one-label file in the same directory and a separate owner's
decision; so is whether the FlexRV annotation group takes `spliceai2_gt040`.

**Four-way at max_MAF 0.001 (2026-10-15, all nine traits; three runs, four arms, same
damaging-missense base, synonymous controls identical across runs -- 5.3 % below 0.05, same
quantiles -- so the differences are the masks).** Gene-trait hits past 2.5e-6 on the
damaging + splice mask: v1@0.2 34, v1@0.5 37, v2@0.25 43, v2@0.40 42; distinct genes across
traits 31 / 34 / 39 / 37, 29 in every arm, only-in-one-arm 1 / 1 / 2 / 2; splice-only hits 0
in every arm except height (0 / 1 / 3 / 1) and alt (0 / 1 / 1 / 0); pooled tail at 1e-3:
3.0e-4 / 2.8e-4 / 2.3e-4 / 2.2e-4, at 1e-4: 1.2e-10 / 6.5e-11 / 4.2e-11 / 2.1e-11 (the
1e-5 quantile at n = 163k is one height gene and is not read). (These are the numbers AFTER
the alt-b5 bin's chr8 was recovered on IMG15 -- the first reading had v2@0.40 at 41 / 36 /
5.9e-10 with GPT's chromosome missing; GPT itself reads 1e-128 on damaging missense alone and
1e-146 with the SpliceAI2 label, 18 decades from the splice calls.) Height carries most of it
(17 / 19 / 22 / 21). Reading: v2 beats v1 in the same direction on every count, and the
ordering v1@0.2 < v1@0.5 reproduces section 16 on independent grounds; the margin
v2@0.40 over v1@0.5 is 4 gene-trait hits and 2 distinct genes, small on correlated traits
(the lipid four share APOB, APOC3, LDLR, ABCA1, LCAT; the anthropometric three share
GIGYF1) but the third contrast in a row to land the same way. Within v2, 0.25 vs 0.40 is a
draw: the band's own cells are null (`v2only@.25` on the synonymous line at every
quantile, 118k cells), the per-trait moves are one gene both ways (0.25 gains height,
bmi, hba1c, alt; loses tchol, t2diab), and the tighter arm has the smaller p in a majority
for every trait. 0.40 stands. **Two metrics retired on this reading:** "best arm by min p
over shared genes" counts mask size (the superset arm v1@0.2 wins every trait by ~1,000
because every gene its band touches and nudges down is a strict win, while a gene no
version touches is a four-way tie credited to all); and the 1e-5 quantile of a 160k-cell
pool is a single gene. The cells: the four-way table and the threshold-count cell are in
`brava-reannotation/workbench/sai2v040_cells.md`'s successor, or reconstructable from the
section-17 reading cell with `MAF = "0.001"` applied inside `table()` and `minp_by_gene()`.

**Touched-gene lift (the last instrument, same day).** For each arm, the gene-traits whose
damaging + splice p differs from the plain `damaging_missense` p in the same run (the splice
set added at least one variant), the combined p against the damaging-only p of the same
genes: touched 145,844 / 129,167 / 136,174 / 121,646 (v1@0.2 / v1@0.5 / v2@0.25 / v2@0.40);
median log10 lift -0.011 / -0.005 / -0.005 / -0.004 and the fraction sharpened 0.48 in every
arm -- for the typical gene every splice set is noise, slightly more dilution than help;
cells below 1e-4 among the touched, combined minus damaging-only: +6 / +8 / +5 / +15. Tighter
beats looser within each version and v2@0.40, touching the fewest genes, gains the most;
counts near 90 carry a sd of ~10, so +15 against +5 is suggestive, but it is the fourth
instrument to order the arms the same way. **Status after this reading (corrected the
same day): v2 over v1 is DECIDED on evidence; 0.25 against 0.40 is a DRAW the data cannot
break** -- counts favour 0.25 (hits 43 vs 41, distinct genes 39 vs 36, the top of the
touched-gene QQ), per-gene behaviour favours 0.40 (the band's own cells null, the typical
gene diluted, lift per touched gene +15 vs +5). The earlier "FINAL 0.40" was the pre-stated
tie-break convention (the ACMG-calibrated cut-off) speaking, not the data, and it was
written as if it were the data. The competing convention is Illumina's README, which names
0.25 as the operating point equivalent to v1's 0.5 -- comparability with every v1-based
result, at no measurable cost on these counts, and a better fit for an inclusion mask than
a clinical evidence-strength threshold. The owner's call, pending one last check: the genes
that are hits under one cut-off and not the other, by name (a cell in the session; few
enough to inspect). If 0.25 is chosen the `out_spliceai2` release files (built at 0.40,
label `spliceai2_gt040`) are rebuilt at 0.25 before upload.

**DECIDED (owner, 2026-10-15): SpliceAI2 at 0.40, on the calibration framing.** The
reannotation repo's own Pejaver-style ClinVar calibration (`brava-reannotation/docs/spliceai.md`,
one rule applied to v1 and SpliceAI2 in the same domain with the same prior, with a simulation,
a REVEL positive control and a permuted negative control) puts v1's Moderate at 0.50 -- the
threshold BRaVa already uses, confirmed by section 16 -- and SpliceAI2's at 0.40; the README's
0.25 is Supporting (0.24) under the same rule, i.e. the README equivalence is recall-matched,
not evidence-matched. So the mask's three predictors (REVEL 0.773, CADD 28.1, SpliceAI2 0.40)
sit at one evidence level by one calibration, which is the citable argument. Soft spots, named:
the number moves with the window and the ClinVar release (0.33-0.50 for SpliceAI2, 0.50-0.66 for
v1), the REVEL control runs one level lenient (circularity), SpliceAI2's training-set overlap
with ClinVar is unchecked; and nine eur traits cannot separate 0.25 from 0.40 on power. The
production mask set keeps a with/without pair -- `pLoF:damaging_missense` against
`pLoF:damaging_missense:spliceai2_gt040` -- so the splice label's contribution is readable in
the production output itself. Release files: the
`out_spliceai2` `new_rules` build (label `spliceai2_gt040`), nothing to rebuild.

## 18. `LOG10P` (2026-10-15): log10 p-value columns in the AoU runs, ON by default

`--is_output_log10P` (`docs/state/handoff-2026-10-14-log10p.md`) adds `log10_Pvalue*`
columns after `SE_Burden` on both region paths and combines the pooled row in log
space, so GPT-class genes (the trait's own enzyme) are estimated rather than floored at 0. The wrapper
(drop-in `harness/work/aou/saige_step2_brava_main_original.sh`) takes `LOG10P`,
**default TRUE** (owner, 2026-10-15); `env["LOG10P"] = "FALSE"` in a submit cell
turns it off, which is REQUIRED for a partial rerun meant for the splice control:
the header differs from the v4 file's and the control refuses. The AoU tools read
columns by name and are unaffected. **Needs IMG15 or later: IMG14 LOST A CHROMOSOME with
the flag on** (GPT for ALT, p ~1e-146 -- a -Inf log twin reached `CCT_log10`, whose error
killed the chromosome from inside `region_postprocess`; `docs/state/handoff-2026-10-14-log10p.md`,
addendum). Any IMG14 job that ran with `LOG10P=TRUE` has its `chrStatus` read and a failed
chromosome rerun on IMG15.

```python
IMG14 = "gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-524a299f"   # --is_output_log10P -- DO NOT USE with LOG10P=TRUE (lost a chromosome, see above)
IMG15 = "gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-f4a316d5"   # the fix; same bytes on Docker Hub as astheeggeggs/saige-slim:1.5.2-dev-f4a316d5 (digest a77c961e); P1 byte-identical, P79 1e-13, the extreme ladder finite in the container
```

The wrapper patch (same shape as section 15's; the three sites):

```bash
# inputs block
LOG10P="${LOG10P:-TRUE}"                              # TRUE | FALSE: --is_output_log10P (group tests; FALSE for a partial rerun to be spliced)

# the command's group block, after the cmd array is built
  if [[ ${LOG10P} == "TRUE" ]]; then cmd+=(--is_output_log10P=TRUE); fi

# the settings echo
LOG10P            = ${LOG10P}
```

## 19. Production rerun on IMG15: the seven-mask SKAT-O, FlexRV, `am_moderate` (2026-10-15)

Decided inputs (owner, 2026-10-15; sections 16-17b): release group files are the `out_spliceai2`
build in `~/Repositories/brava-reannotation/out_spliceai2/aou/` (SpliceAI2 at 0.40, label
`spliceai2_gt040`); 23 chromosomes each of `new_rules`, `flexrv_AM`, `am_moderate` (chrX is
uploaded and not run). Image IMG15. `LOG10P` on (the wrapper default). Masks:

```python
IMG15 = "gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-f4a316d5"
REL = f"{BUCKET}/SAIGE_SLIM_TESTING/annotations_v5"          # the release upload directory
PROD_MASKS = ("pLoF,damaging_missense,pLoF:damaging_missense,pLoF:damaging_missense:spliceai2_gt040,"
              "synonymous,other_missense,pLoF:damaging_missense:spliceai2_gt040:other_missense:synonymous")
FLEXRV_GROUP = "pLoF:damaging_missense:other_missense:spliceai2_gt040"   # the splice bin holds missense with a splice call; AM-scored
AM_MASK = "am_moderate"                                                 # AlphaMissense >= 0.906, missense only (checked: every variant scores >= 0.906 in the FlexRV file, the only scored >= 0.906 outside it are pLoF)
nr_grp = lambda c: f"{REL}/aou.chr{c}.new_rules.saige_group.txt.gz"
fx_grp = lambda c: f"{REL}/aou.chr{c}.flexrv_AM.txt.gz"
am_grp5 = lambda c: f"{REL}/aou.chr{c}.am_moderate.saige_group.txt.gz"
```

### 19.1 Local prep (laptop): one directory, chr1 `new_rules` split, md5 of everything

Only `aou.chr1.new_rules.saige_group.txt.gz` (19.0 MB) is above the size that uploaded intact in
section 17 (15 MiB parts); everything else goes up whole.

```bash
cd ~/Repositories/brava-reannotation/out_spliceai2/aou
mkdir -p release_v5 && cp aou.chr*.new_rules.saige_group.txt.gz aou.chr*.flexrv_AM.txt.gz aou.chr*.am_moderate.saige_group.txt.gz release_v5/
cd release_v5 && split -b 15m aou.chr1.new_rules.saige_group.txt.gz aou.chr1.new_rules.saige_group.txt.gz.part_ && rm aou.chr1.new_rules.saige_group.txt.gz
ls -l aou.chr1.new_rules.saige_group.txt.gz.part_*          # two parts, 15,728,640 + the rest
cd .. && (cd release_v5 && for f in *; do printf '%s %s %s\n' "$(md5 -q "$f")" "$(stat -f %z "$f")" "$f"; done) > release_v5.md5.txt
printf '%s %s %s\n' "$(md5 -q aou.chr1.new_rules.saige_group.txt.gz)" "$(stat -f %z aou.chr1.new_rules.saige_group.txt.gz)" aou.chr1.new_rules.saige_group.txt.gz >> release_v5.md5.txt
wc -l release_v5.md5.txt                                      # 71: 68 whole files + 2 parts + the composed chr1
```

Upload `release_v5/` (70 objects) and `release_v5.md5.txt` into `SAIGE_SLIM_TESTING/annotations_v5/`.

### 19.2 Workbench: compose chr1, md5-check all 69 objects

```python
import hashlib
D = REL
parts = [f"{D}/aou.chr1.new_rules.saige_group.txt.gz.part_{x}" for x in ("aa", "ab")]
dest = f"{D}/aou.chr1.new_rules.saige_group.txt.gz"
sizes = {}
for p in parts:
    out = subprocess.run(["gcloud", "storage", "ls", "-l", p], capture_output=True, text=True).stdout.split()
    sizes[p] = int(out[0]) if out else -1
print(sizes)                                                   # part_aa must be 15728640 exactly
assert sizes[parts[0]] == 15728640, "part_aa truncated: re-upload it"
print(subprocess.run(["gcloud", "storage", "objects", "compose", *parts, dest], capture_output=True, text=True).stderr[-200:])
ref = {}
for l in subprocess.run(["gcloud", "storage", "cat", f"{D}/release_v5.md5.txt"], capture_output=True, text=True).stdout.splitlines():
    h, n, f = l.split(); ref[f] = (h, int(n))
bad = []
for f, (h, n) in sorted(ref.items()):
    if ".part_" in f: continue
    data = subprocess.run(["gcloud", "storage", "cat", f"{D}/{f}"], capture_output=True).stdout
    ok = hashlib.md5(data).hexdigest() == h and len(data) == n
    if not ok: bad.append(f)
    print(("OK  " if ok else "BAD ") + f, len(data))
print(len(ref) - 2, "objects checked;", "ALL OK" if not bad else f"RE-UPLOAD: {bad}")
```

Nothing below runs until that prints `ALL OK`.

### 19.3 Two canaries on IMG15 (on demand), read before the fleet

One eur quantitative SKAT-O bin and the matching FlexRV trait. What to read in each: `chrStatus`
`failed none`; the log's settings echo `LOG10P            = TRUE` and the `Running` line carrying
`--is_output_log10P=TRUE`; the results header with `log10_Pvalue` after `SE_Burden`; the end-of-run
lines (`SPA fallbacks`, `Exact test by weight`, `Stretch gate`, `SKAT-O cells` / `BURDEN cells`);
and for FlexRV the one-line `.singleAssoc.txt` placeholder (`SINGLE_IN_GROUP=FALSE`).

```python
load_runs()
v4 = next(k for k in bins_of("eur", "ldlc") if 19 in RUNS[k]["chrs"])        # the LDLR bin
r = RUNS[v4]; env = {k: v for k, v in r["env"].items() if k != "REGIONS"}
submit(name="canary-sk6-eur-ldlc", image=IMG15, masks=PROD_MASKS, fasttest=r["fasttest"], outdir="eur/canary_v6",
       anc="eur", pheno=r["pheno"], trait=r["trait"], firth=r["firth"], env=env, geno=r["geno"], chrs=r["chrs"],
       plk=lambda c, a="eur": cod(c, a), grp=nr_grp, testtype="group", machine_type=r["machine_type"],
       spot=False, pooled=True, sidecars=True, sparse_grm=r["sparse_grm"])
f4 = next(k for k in RUNS if k.startswith("fxv4-eur-ldlc"))
r = RUNS[f4]; env = {k: v for k, v in r["env"].items() if k != "REGIONS"}; env["SINGLE_IN_GROUP"] = "FALSE"
submit(name="canary-fx6-eur-ldlc", image=IMG15, masks=FLEXRV_GROUP, fasttest=r["fasttest"], outdir="eur/canary_v6",
       anc="eur", pheno=r["pheno"], trait=r["trait"], firth=r["firth"], env=env, geno=r["geno"], chrs=r["chrs"],
       plk=lambda c, a="eur": cod(c, a), grp=fx_grp, testtype="group", machine_type=r["machine_type"],
       spot=False, pooled=True, sidecars=True, sparse_grm=r["sparse_grm"])
```

```python
# reading the canaries
for nm in ("canary-sk6-eur-ldlc", "canary-fx6-eur-ldlc"):
    pre = RUNS[nm]["out_prefix"]; log = fetch_log(nm)
    print("==", nm)
    print("  chrStatus:", rows(pre + ".chrStatus.txt.gz"))
    print("  LOG10P echo:", [l for l in log.splitlines() if l.startswith("LOG10P")][:1], "| flag in cmd:", "--is_output_log10P=TRUE" in log)
    h = rows(pre + ".txt.gz")[0].split("\t"); print("  log10 columns:", [c for c in h if c.startswith("log10_")], "after SE_Burden:", h[h.index("SE_Burden") + 1] if "SE_Burden" in h else None)
    for key in ("SPA fallbacks", "Exact test by weight", "Stretch gate", "SKAT-O cells", "BURDEN cells", "Concentrated cells"):
        for l in log.splitlines():
            if l.startswith(key): print("  ", l[:160]); break
```

### 19.4 The fleet: every v4 record, same bins, new image and files

Three job sets, each derived from the v4 registry so bins, models, `fasttest`, `firth`, `geno`,
machine type and the GRM choice are the v4 job's; only name, image, masks, group file and output
directory change (plus `SINGLE_IN_GROUP=FALSE` on FlexRV AND on `am_moderate`: the single-variant
file is written once, by the seven-mask run, whose variants are a superset). **Launched 2026-10-15
(owner): Spot for every job** (the standing preference; section 16 lost 2 of 80, the mop-up pass
recovers losses at the cost of their wasted hours) **and the v4 machine types: `n2-highmem-2`
for eur, `n2-standard-2` for the rest.** (A standard-2 for eur was proposed and withdrawn the same
hour: the 2026-09-25 sizing handoff has eur SKAT-O at 6,057-6,170 MB, flat across chromosomes --
a property of the model, not the work -- which is why eur is on highmem; the chr1 1,188 / chr2
2,134 / chr19 2,230 MB figures in the wrapper are AFR's. Today's eur canary reported a SAMPLED
`PEAK_RSS_MB` of 1,191 on IMG15, five times below September's; unexplained, not relied on.)
Canaries read clean first: `chrStatus` none failed, `LOG10P` echoed and in the command, the
log10 columns after `SE_Burden`, the single-variant file real on sk6 and the placeholder on fx6/am6,
LDLR pLoF 4.4e-10 / damaging 9.5e-19 (v4's) / pLoF:damaging 3.1e-26 / with the splice label 2.2e-26
/ synonymous 1 / am_moderate 9.1e-16; SAIGE seconds 1,012 (sk6, two eur chromosomes), 165 (fx6, eur
chr1, 315k BURDEN cells), 110 (am6). 806 jobs per set, 2,418 in all; one set armed per cell
(submission is ~1 h per set), `save_runs()` after each. `ARMED = False` prints the plan.

```python
ARMED = False
def fleet(kind):
    """kind: 'sk6' (SKAT-O, PROD_MASKS, new_rules), 'fx6' (FlexRV, FLEXRV_GROUP, flexrv_AM), 'am6' (am_moderate)."""
    src = "fxv4-" if kind == "fx6" else "skv4-"
    plan = []
    for v4 in sorted(k for k in RUNS if k.startswith(src)):
        r = RUNS[v4]; anc = r["anc"]
        if kind == "sk6":  masks, grp, od = PROD_MASKS, nr_grp, f"{anc}/skato_v6"
        if kind == "fx6":  masks, grp, od = FLEXRV_GROUP, fx_grp, f"{anc}/flexrv_v6"
        if kind == "am6":  masks, grp, od = AM_MASK, am_grp5, f"{anc}/am_v6"
        nm = v4.replace(src, kind + "-")
        if already(out_prefix_for(nm, anc, r["pheno"], r["chrs"]), od): continue
        env = {k: v for k, v in r["env"].items() if k != "REGIONS"}
        if kind in ("fx6", "am6"): env["SINGLE_IN_GROUP"] = "FALSE"   # every am_moderate variant is a new_rules missense variant: the sk6 run's single-variant file covers both (owner, 2026-10-15)
        spot = True                          # the mop-up pass sets False for the stragglers
        machine_type = "n2-highmem-2" if anc == "eur" else "n2-standard-2"   # the v4 sizing (see above)
        plan.append((nm, anc, r["pheno"], r["chrs"], spot))
        if ARMED:
            submit(name=nm, image=IMG15, masks=masks, fasttest=r["fasttest"], outdir=od, anc=anc, pheno=r["pheno"],
                   trait=r["trait"], firth=r["firth"], env=env, geno=r["geno"], chrs=r["chrs"],
                   plk=lambda c, a=anc: cod(c, a), grp=grp, testtype="group", machine_type=machine_type,
                   spot=spot, pooled=True, sidecars=True, sparse_grm=r["sparse_grm"])
    print(kind, len(plan), "jobs", "(submitted)" if ARMED else "(plan)", "| spot", sum(p[4] for p in plan), "on demand", sum(not p[4] for p in plan))
    return plan
plan_sk = fleet("sk6"); plan_fx = fleet("fx6"); plan_am = fleet("am6")
if ARMED: save_runs()
```

Run unarmed first: the three counts must equal the v4 registry's (every skv4 record twice, every
fxv4 once). Arm, run, `save_runs()`. The same cell IS the mop-up: `already()` keys on the delivered
results file, so rerunning it after the batch settles resubmits exactly the Spot losses and
nothing else (set `spot = False` inside for that pass).

### 19.4b The first batch's failures (2026-10-15 evening) and the mop-up

35 of 2,418 failed, none labelled a preemption: Batch's `--full` events on this provider carry
only `scheduled` / `start` / `fail` and the status message is the word FAILURE, so a filter on
"preempt" or exit code 50001 finds nothing. The signature IS the preemption: empty stderr,
stdout cut off mid-run (one at `T_SAIGE_START`, two in the middle of their output, four amr
FlexRV jobs with nothing written at all), no SAIGE error anywhere. The failures cluster by trait
(all ten eur ifhern bins, all four amr t2diab FlexRV bins) because a trait's bins are submitted
within seconds, land in one zone and are reclaimed together -- the zone correlation. One
(`sk6-eur-height-b5`) had finished SAIGE and was killed in delivery. `fetch_log` reads the
runner log (`<job>.log` sorts after `-stderr.log`); the SAIGE output is the `-stdout.log`.
**Mop-up**: wait until nothing is RUNNING, then `_DELIVERED.clear()` (the `delivered()` cache
was taken when nothing had been delivered -- without clearing it the pass resubmits all 806),
`spot = False` in `fleet()`, unarmed to see the list, armed per set, `save_runs()`.

### 19.4c The summary job's chromosome scan found four lost chromosomes (2026-10-16): LEDGER #181

afr HTN (`sk6-afr-htn-b2` chr 5, 20; `am6-afr-htn-b0/b1/b2` chr 6, 12, 5) and eur HTN (`sk6-eur-htn-b3` chr 11); `Error in if
(VarQ_2 == 0)` in the region path's burden rescale, a NaN from the R-side saddlepoint outside
the per-cell catch (`docs/state/handoff-2026-10-16-region-spa-na.md`). Fixed the same day
(the NaN is dumped, counted, named, and the rescale takes the normal p as `SPA_pval` itself
does on non-convergence); image IMG16. The four jobs are rerun in full on IMG16 with
`SPA_DUMP=1` (the dumps say which kernel path made the NaN), then afr's summary once more.
Every other ancestry's scan is read the same way before its pass 2. **Was anything else hit?**
The kernel's NaN on the C++ test paths would have printed an `NA` p rather than crashing; the
check is per (trait, gene) on a mask whose definition is the same in v4 and v6 (`pLoF`): cells
non-`NA` in v4 and `NA` in v6, binary traits against the quantitative baseline (which has no
saddlepoint). sas: 0 such cells in every trait (2026-10-16). (A per-trait NA COUNT comparison
is confounded by the annotation rebuild -- genes gain and lose cells -- and a raw NA fraction
by the mask composition; both were tried first and read as alarming before the per-cell
statistic was written.)

```python
IMG16 = "gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-3a6d81fa"   # LEDGER #181; same bytes on Docker Hub (digest 4b987d8f); P1 byte-identical and the planted NA net checked in the container
redo = ["sk6-afr-htn-b2", "am6-afr-htn-b0", "am6-afr-htn-b1", "am6-afr-htn-b2", "sk6-eur-htn-b3"]
for nm in redo:
    r = RUNS[nm]; env = dict(r["env"]); env["SPA_DUMP"] = "1"
    kind = nm[:3]; grp = {"sk6": nr_grp, "am6": am_grp5, "fx6": fx_grp}[kind]
    submit(name=nm, image=IMG16, masks=r["masks"], fasttest=r["fasttest"], outdir=r["outdir"], anc=r["anc"],
           pheno=r["pheno"], trait=r["trait"], firth=r["firth"], env=env, geno=r["geno"], chrs=r["chrs"],
           plk=lambda c, a=r["anc"]: cod(c, a), grp=grp, testtype="group", machine_type=r["machine_type"],
           spot=False, pooled=True, sidecars=True, sparse_grm=r["sparse_grm"])
save_runs()
```

**The sidecar's delivery (added after the reruns were submitted, so THEY carry only the log
line and the dumps):** the wrapper drop-in now takes `OUT_FILE_SPAREGIONNA` (delivered like the
other sidecars, a one-line placeholder when the run had none), and `submit()`'s sidecar block
needs one more line beside `OUT_FILE_STRETCHGATE`:
`"--output", f"OUT_FILE_SPAREGIONNA={od}/{pre}.spaRegionNA.txt.gz",`. The summary tool scans
the sidecar on all three arms, prints `spaRegionNA: N rescale saddlepoint(s) returned nothing
...; <arm>/<trait> K` and writes `<anc>.spaRegionNA.<tag>.tsv.gz` (selftest 31 checks).
**IMG17 (`1.5.2-dev-b5456dee`, same day): the counter was blind on IMG16** -- the kernel guard
caught the NaN first and R cannot see `SPA_pval`'s convergence flag (passed by value), so the
five reruns read `0` and no dump fired. `SPA_pval_status` fixes that: every non-converged
rescale is now counted with its reason and dumped. The cause is read from ONE bin rerun on
IMG17 with the dump on (the production outputs from IMG16 are complete and correct; this rerun
is the instrument, its output overwrites the same prefix with identical cells):

```python
IMG17 = "gcr.io/weighty-elf-452116-c7/saige-slim:1.5.2-dev-b5456dee"
nm = "sk6-afr-htn-b2"; r = RUNS[nm]; env = dict(r["env"]); env["SPA_DUMP"] = "1"
submit(name=nm, image=IMG17, masks=r["masks"], fasttest=r["fasttest"], outdir=r["outdir"], anc=r["anc"],
       pheno=r["pheno"], trait=r["trait"], firth=r["firth"], env=env, geno=r["geno"], chrs=r["chrs"],
       plk=lambda c, a=r["anc"]: cod(c, a), grp=nr_grp, testtype="group", machine_type=r["machine_type"],
       spot=False, pooled=True, sidecars=True, sparse_grm=r["sparse_grm"])
```

Then its log's `Region saddlepoint NA: N ... K returned a non-finite value` line and, in the
delivered `spaDump.tgz`, the `spaRegionNA_*.tsv` files (header: q, qinv, the normal log p,
n_dropped, n; columns mu, g). Read them where the data lives: a vector whose non-zero `g` are
all cases with mu near 1, or a root at the attainable edge, is the #168-branch shape.
Read each rerun's log for `Region saddlepoint NA: N ...` and `chrStatus` `failed none`; then
`aou_spa_dumpcheck.sh` is NOT the tool for the R dumps (different header) -- read them where
the data lives: `n_dropped` in the header says whether the pinned guard emptied the vector,
and a vector with no non-zero `g` says the same.

### 19.5 Before any reading: the chromosome scan -- NOW DONE BY PASS 1 (2026-10-16)

The notebook version below killed a session (2,418 sidecars at ~2 s each through `gcloud cat`,
the browser connection dropped and the output with it). `aou_flexrv_summary.sh` now scans every
job's `chrStatus` on the VM, writes `<anc>.chrstatus.v6.tsv.gz` and prints ONE line in the job
log -- `chrStatus: N jobs scanned, K with a failed chromosome[: arm/job(chr c) ...]` -- which is
read before that ancestry's summaries are used. The notebook cell is kept for a single
ancestry or a spot check only.

```python
for od_suffix in ("skato_v6", "flexrv_v6", "am_v6"):
    for anc in ("afr", "amr", "eas", "eur", "sas"):
        ls_ = subprocess.run(["gcloud", "storage", "ls", f"{BUCKET}/SAIGE_SLIM_TESTING/{anc}/{od_suffix}/*.chrStatus.txt.gz"], capture_output=True, text=True).stdout.split()
        bad = [(p.rsplit("/", 1)[-1].split(".")[0], dict(l.split("\t") for l in rows(p)).get("failed")) for p in ls_]
        bad = [b for b in bad if b[1] not in (None, "none")]
        print(f"{anc}/{od_suffix}: {len(ls_)} jobs, {len(bad)} with a failed chromosome", bad[:5])
```

A failed chromosome is a job to rerun in full (same name, same bins) before its trait is read;
both lost-chromosome episodes (#170, the IMG14 log10 defect) passed every other check.

### 19.6 Reading (morning of 2026-10-16): three dsub passes, nothing read into the notebook

The reading follows the v4 pattern (section 7): gene-level summaries are built INSIDE dsub
jobs, figures are rendered INSIDE dsub jobs, and the notebook only displays PNGs and a text
file. Every number is a count of genes or cells. The three tools, each with a fabricated-data
selftest run before this was written (`bash harness/tools/<tool>_selftest.sh` locally; the
python ones inside the pyplot image: `docker run --rm -v $PWD:/repo -w /repo pyplot:20260928
bash harness/tools/aou_v6_vs_v4_selftest.sh python3`):

- `aou_flexrv_summary.sh` (RUNTAG=v6): reads `skato_v6`, `flexrv_v6` AND `am_v6` (the
  am_moderate cells go into the wide cells table as `am_moderate@<maf>`, never into the
  SKAT-O per-gene best), writes `<anc>.flexrv_vs_skato.v6.tsv.gz`, `<anc>.skato_cells.v6.tsv.gz`
  and NEW `<anc>.floor.v6.tsv.gz` -- every cell past log10 p = -300 by its log10 column, the
  genes the flag exists for (empty for v4). Selftest 21 checks.
- `aou_flexrv_qq.py`: the per-(ancestry, trait) QQ grids, now over whatever masks the cells
  table holds (v6: the seven plus am_moderate, a 3 x 4 grid). Selftest 6 checks.
- `aou_v6_vs_v4.py` (NEW): one ancestry, the two summary sets; prints and writes
  `<anc>.v6_vs_v4.txt` with the CONTROLS first (synonymous below 0.05 old vs new; pLoF gene
  for gene as the cross-image control: Spearman, cells moved by more than a decade -- expected
  only past 1e-5 -- and the hit counts), then the contrasts (3. the splice label; 4a/4b the
  reannotation like for like and missense-only; 5. am_moderate against damaging_missense;
  6. FlexRV pooled against the splice mask's SKAT-O p), a per-trait table
  `<anc>.v6_vs_v4.tsv`, and two figures `<anc>.controls.png`, `<anc>.contrasts.png`.
  Selftest 11 checks (planted moves, splice hits, am hits, FlexRV hits all recovered).

Copy the three tools next to `aou_dsub.bash` on the workbench first (they are `--script`
arguments). The v4 summaries already exist under `SAIGE_SLIM_TESTING/summary_v4/` (section 7).

```python
ANCS = ["afr", "amr", "eas", "eur", "sas"]
SV6 = f"{BUCKET}/SAIGE_SLIM_TESTING/summary_v6"; SV4 = f"{BUCKET}/SAIGE_SLIM_TESTING/summary_v4"
PYPLOT = "gcr.io/weighty-elf-452116-c7/pyplot:20260928"
def run_dsub(args):
    script = "source /home/jupyter/aou_dsub.bash\n" + " ".join(shlex.quote(a) for a in args)
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=os.environ)
    print(args[2], "rc =", out.returncode, "|", ((out.stdout or "") + (out.stderr or "")).strip()[-300:])
# pass 1: the v6 summaries, one job per ancestry (eur ~1 h on the cloud-cli image: 806 x 3 files streamed)
for anc in ANCS:
    run_dsub(["aou_dsub", "--name", f"summary6-{anc}", "--image", "gcr.io/google.com/cloudsdktool/google-cloud-cli:slim",
              "--env", f"ANC={anc}", "--env", f"BASE={BUCKET}/SAIGE_SLIM_TESTING", "--env", "RUNTAG=v6",
              "--output-recursive", f"OUT_DIR={SV6}", "--machine-type", "n2-standard-2", "--disk-size", "20",
              "--script", "aou_flexrv_summary.sh"])
```

Read each summary job's log before the next pass: the lines `flexrv: N v4 files ...; skato: M
files`, `am_moderate: K files` (806 / 806 / 806 for eur... per ancestry the registry's count)
and `floor table: F row(s)`.

```python
# pass 2 (after pass 1 delivers): QQ grids and the v6-vs-v4 comparison, per ancestry
for anc in ANCS:
    run_dsub(["aou_dsub", "--name", f"qq6-{anc}", "--image", PYPLOT, "--env", f"ANC={anc}",
              "--input", f"IN_SUMMARY={SV6}/{anc}.flexrv_vs_skato.v6.tsv.gz", "--input", f"IN_CELLS={SV6}/{anc}.skato_cells.v6.tsv.gz",
              "--output-recursive", f"OUT_DIR={SV6}/qq/{anc}", "--machine-type", "n2-standard-2", "--disk-size", "20",
              "--script", "aou_flexrv_qq.py"])
    run_dsub(["aou_dsub", "--name", f"cmp6-{anc}", "--image", PYPLOT, "--env", f"ANC={anc}",
              "--input", f"IN_SUM_NEW={SV6}/{anc}.flexrv_vs_skato.v6.tsv.gz", "--input", f"IN_CELLS_NEW={SV6}/{anc}.skato_cells.v6.tsv.gz",
              "--input", f"IN_FLOOR_NEW={SV6}/{anc}.floor.v6.tsv.gz",
              "--input", f"IN_SUM_OLD={SV4}/{anc}.flexrv_vs_skato.v4.tsv.gz", "--input", f"IN_CELLS_OLD={SV4}/{anc}.skato_cells.v4.tsv.gz",
              "--output-recursive", f"OUT_DIR={SV6}/cmp/{anc}", "--machine-type", "n2-standard-2", "--disk-size", "20",
              "--script", "aou_v6_vs_v4.py"])
```

```python
# reading: the text report first (controls at the top), then the two figures per ancestry
from IPython.display import Image, display
for anc in ANCS:
    print(subprocess.run(["gcloud", "storage", "cat", f"{SV6}/cmp/{anc}/{anc}.v6_vs_v4.txt"], capture_output=True, text=True).stdout[:6000])
    for fig in ("controls", "contrasts"):
        subprocess.run(["gcloud", "storage", "cp", f"{SV6}/cmp/{anc}/{anc}.{fig}.png", f"/home/jupyter/{anc}.{fig}.png"], check=True)
        display(Image(f"/home/jupyter/{anc}.{fig}.png"))
```

How to read it, decided before the numbers:
1. `synonymous below 0.05` must be ~0.05 in both runs and the synonymous QQ on the band;
   if the new run inflates where the old did not, stop -- that is the image or the model.
2. pLoF: Spearman near 1; the cells moved by more than a decade should be few and past
   1e-5 (the saddlepoint tail and exact-by-weight changes between the images); hit counts
   old vs new within a handful. Many moved cells above 1e-5 means the arithmetic changed
   where it must not have, and that is read before any contrast.
3. The splice label (3.): the only-second count is what SpliceAI2 at 0.40 adds; only-first
   is what it dilutes away. From the eur evidence run the expectation is a small positive net
   with the band of borderline swaps both ways.
4. The reannotation (4a/4b): the net change of the damaging mask; distinct genes per ancestry.
5. am_moderate (5.) and FlexRV (6.): what each buys over the damaging mask, as hit counts and
   the scatter's spread above the identity line.
6. The floor section: the genes past 1e-300, by name and log10 -- the ones a linear p could
   not rank; the known shape is one or a few cis-enzyme genes per trait (GPT for ALT).
Then the per-(ancestry, trait) QQ grids under `summary_v6/qq/<anc>/` for anything the
report flags, exactly as in section 7.

### 19.7 The reading (2026-10-16, afternoon): controls clean on five ancestries, the gains by size

Every number below is a count of genes or cells from `<anc>.v6_vs_v4.txt` (cells at max_MAF
0.001, hits past 2.5e-6, gene for gene on the shared (trait, gene) pairs). The five reports and
the figures are under `summary_v6/cmp/<anc>/`, the per-(ancestry, trait) QQ grids under
`summary_v6/qq/<anc>/`.

**Controls (read first; all five clean).**
- synonymous below 0.05: sas 0.048 / 0.048, eas 0.048 / 0.048, amr 0.051 / 0.051, afr 0.051 /
  0.051, eur 0.051 / 0.051 (old / new, identical to three digits; eur on the same 791,559 cells).
  Per trait, the polygenic ones sit where they sat in v4 (amr Height 0.075, BMI 0.064).
- pLoF gene for gene (the cross-image control, the one mask whose definition did not change):
  Spearman 1.0000 in all five; cells moved by more than a decade 0 / 0 / 0 / 1 / 0; hits and
  distinct genes identical (eur 85 hits, 67 genes; amr and afr 27 / 23; eas 5; sas 1). The two
  images agree on the unchanged mask exactly: everything below is annotation and method.
- The chromosome scans: 0 failed everywhere after the HTN reruns (19.4c). The `NA` check per
  (trait, gene) on pLoF: no cell v4 tested came back `NA` in v6 outside HTN's lost chromosomes.

**The contrasts, by size.**
- Reannotation, like for like (old `pLoF;damaging_missense` with SpliceAI v1 >= 0.2 inside, vs
  new `pLoF;damaging_missense;spliceai2_gt040`): a net gain wherever there is power -- amr 22 ->
  26 hits (19 -> 23 genes), afr 27 -> 31 (23 -> 26), eur 95 -> 110 (80 -> 89 genes; 26 gained,
  11 lost); sas 4 -> 5, eas 4 -> 3. Spearman ~0.5: a different mask (REVEL 0.773 / CADD 28.1 /
  the AlphaMissense fill), tens of thousands of cells moved by a decade.
- The splice label (new `pLoF;damaging_missense` vs the same with `spliceai2_gt040`): marginal
  below eur (afr 32 vs 31, amr 26 vs 27, eas 4 vs 3, sas 6 vs 6), a real net in eur: 102 -> 110
  hits, 84 -> 89 distinct genes (13 gained, 5 lost). And on the one gene where splice calls
  matter most it is 19 decades: GPT for ALT, log10 p -383 with the label, -365 without.
- `am_moderate` vs `damaging_missense`: fewer hits alone (eur 34 vs 44, afr 10 vs 13, amr 10 vs
  13) but a partly different set (eur 14 only it finds): complementary, not redundant.
- FlexRV pooled vs the splice mask's SKAT-O p: the largest margin in every cohort with power --
  eur 138 vs 110 hits (49 only FlexRV, 21 only SKAT-O; 110 vs 89 distinct genes; Height 52 vs
  35, HDLC 18 vs 14, TG 13 vs 10, T2Diab 7 vs 6), afr 41 vs 31 (Height 14 vs 8), amr 41 vs 27
  (Height 12 vs 5), eas 6 vs 3, sas 4 vs 6. **Caveat that stays until it is measured: FlexRV has
  no null control of its own in this run** (no synonymous arm; its calibration rests on the sim
  study and on the per-trait QQ grids, where its pooled panel sits on the band on null traits).
  The control proposed: a `flexrv_RAND` group file (the `score:AM` line replaced by seeded
  uniform values) run on `synonymous` (the machinery's null) and on the damaging group (the
  score's positive control against `AM`), eur, the nine evidence traits.
- The floor (`<anc>.floor.v6.tsv.gz`): nothing past log10 p = -300 in sas, eas, amr, afr; eur has
  four cells, all GPT for ALT at max_MAF 0.01 (-383 splice mask, -365 missense-only, -342 all-in,
  Cauchy -382). A linear `Pvalue` reads 0 for all four; `--is_output_log10P` ranks them.

**What this run cost in corrections**: 35 Spot preemptions resubmitted (19.4b); five HTN jobs
rerun for LEDGER #181 (19.4c); one stray `fxv5` run excluded from every ancestry's FlexRV
directory by the summary tool's prefix guard; the notebook chromosome scan replaced by the
summary job's.

### 19.8 FlexRV outliers (2026-10-16): the sparse-rung shape on quantitative traits

For eur's 138 FlexRV hits, the decades FlexRV sits below SKAT-O's best cell (`gap`), the best
set and whether the pooled p is one set or many. The large gaps are textbook genes sharpened
(ABCA1 +21 decades, LDLR +15, LCAT +14, ANGPTL3, PCSK9, LPL, PKD1, FBN1, IHH, ADAMTS10, BRCA1,
MLH1), all exome-wide in SKAT-O too, best sets spread over the transforms. Of the 21 FlexRV-only
hits, 19 have a SKAT-O cell between 2.5e-6 and 1e-3 (CFH for AMD at 1e-3 -> 1e-7, LDLR for CAD,
IGF1R, ANTXR1, PPARG, NR1H3, INSR among them) and two have no SKAT-O signal at all.

**The shape to hold back**: a hit carried by a sparse rung (`MACsparse` / `MACsparse` / `MACsparse`) with a
steep score weighting on a QUANTITATIVE trait. `sbeta_b` is `s^(b-1)` -- `sbeta20` is zero
below AlphaMissense 0.67 after rounding -- so such a set is its top-scored singletons, and on an
inverse-normalised trait one carrier contributes at most z ~ 4.9 at eur's N: two or three tall
carriers reach 1e-6..1e-8 with no gene-level effect. Checked by Kish's effective carrier count,
(sum w)^2 / sum w^2 over the carriers' weights, from the group file's `score:AM` line and the
`sk6` run's per-variant allele counts (the cell in the session; reported as a BAND, it derives
from real-data carrier counts): three eur Height hits -- ENSG00000186970 (`sbeta20_MACsparse`,
SKAT-O 0.12), ENSG00000140564 (FURIN, `sbeta20_MACsparse`, SKAT-O 1.5e-4), ENSG00000159224
(`cbrt0.8_MACsparse`, SKAT-O 1.3e-3) -- have FEWER THAN FIVE effective carriers each, with 20+
variants and MAC >= 20 on the set (the plain counts do not show it; the weighting does). Held
back. SRSF2 for HDL-C (`sbeta10_mafbeta5`, SKAT-O 1.8e-4) has 20+ effective carriers and stands
as a borderline FlexRV-only hit of unknown biology.

**Proposed (owner's call, not built):** emit Kish's effective carrier count per set as a column on
the BURDEN path for every trait (it is already computed for the binary exact-by-weight routing,
LEDGER #172) and filter FlexRV hits below ~10 at reporting; or drop the `MAC<=1`/`MAC<=3` rungs
on quantitative traits. The random-score null (19.7) sizes how often the rungs do this under no
signal, which is the number that decides between the two.
