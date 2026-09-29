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
    print(subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=os.environ).stdout[-300:])
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
absence of `ERROR`. The notebook's `aou_spa_splice.sh` must be the repo's
current version (self-test: 3 positive controls, 7 planted refusals).

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

If a count is 0 and `STEP1` names another directory, put that directory in
`MODEL_DIR`. Then:

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
