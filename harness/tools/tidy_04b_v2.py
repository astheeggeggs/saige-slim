#!/usr/bin/env python3
"""Build 04b_production_tidy_v2.ipynb from NEW_04b_production_tidy.ipynb.

Three differences from the 2026-09-29 builder:

  1. Cells are selected by a substring UNIQUE ANYWHERE IN THE SOURCE, not by
     the first line. Sixteen first lines collide (six of them function
     definitions), and the four `submit_check_v7` cells are indistinguishable
     by signature -- they differ only in the body, which is where the trap-14a
     bug lives.
  2. `defonly()` copies a cell up to a stop line, so a definition can be kept
     while the one-off driver loop stapled to its bottom is dropped. Three
     `submit_check_v7` cells, `missing_v7`, `stamp_for`, `classify` and six
     others are def+driver pairs.
  3. Two surgical, ASSERTED replacements fix live bugs (see FIXES below). A
     replacement that does not match is a hard error, never a silent no-op.

Nothing is retyped: every kept cell is the notebook's own bytes.
Run in the workbench:   python3 tidy_04b_v2.py
The original notebook is not modified.
"""
import nbformat, sys

OLD = "NEW_04b_production_tidy.ipynb"
NEW = "04b_production_tidy_v2.ipynb"

nb   = nbformat.read(OLD, as_version=4)
code = [c for c in nb.cells if c.cell_type == "code"]
_used, _report = [], []


_seen = {}


def _claim(cell, marker):
    """Refuse to emit one source cell twice. tail_log, numdiff and twenty other
    helpers live INSIDE the 398-line submit() cell, so a marker naming one of
    them silently selects that whole cell a second time."""
    key = id(cell)
    if key in _seen:
        sys.exit(f"CELL ALREADY TAKEN by {_seen[key]!r}; {marker!r} names the same "
                 "cell -- it is probably a helper defined inside a larger one")
    _seen[key] = marker


def _find(marker, allow_dupes=False):
    hits = [c for c in code if marker in (c.source or "")]
    if not hits:
        sys.exit(f"MARKER NOT FOUND: {marker!r}")
    if len(hits) > 1:
        if not allow_dupes:
            idx = [i for i, c in enumerate(code) if marker in (c.source or "")]
            sys.exit(f"MARKER MATCHES {len(hits)} CELLS (code-cell indices {idx}): "
                     f"{marker!r} -- give a longer one")
        srcs = {h.source.rstrip() for h in hits}
        if len(srcs) != 1:
            sys.exit(f"{marker!r} matches {len(hits)} cells that are NOT identical "
                     "-- they must be disambiguated by hand, not deduped")
    return hits[0]


def old(marker, allow_dupes=False):
    """A kept cell, verbatim."""
    c = _find(marker, allow_dupes)
    _claim(c, marker)
    _used.append(marker)
    return nbformat.v4.new_code_cell(c.source.rstrip())


def defonly(marker, stop, allow_dupes=False):
    """A kept cell up to (not including) the first line containing `stop`:
    the definition without the one-off driver loop underneath it."""
    c = _find(marker, allow_dupes)
    lines, cut = c.source.splitlines(), None
    for n, l in enumerate(lines):
        if stop in l:
            cut = n
            break
    if cut is None:
        sys.exit(f"STOP LINE NOT FOUND in the cell matching {marker!r}: {stop!r}")
    if cut == 0:
        sys.exit(f"stop line is the FIRST line of {marker!r} -- nothing would be kept")
    _claim(c, marker)
    _used.append(marker)
    _report.append(f"  defonly {marker[:52]:54s} kept {cut}/{len(lines)} lines")
    return nbformat.v4.new_code_cell("\n".join(lines[:cut]).rstrip())


def patched(marker, subs):
    """A kept cell with asserted literal replacements. subs = [(before, after)]."""
    c = _find(marker)
    src = c.source
    for before, after in subs:
        n = src.count(before)
        if n != 1:
            sys.exit(f"PATCH MATCHED {n} TIMES (need exactly 1) in {marker!r}:\n"
                     f"  {before[:160]!r}")
        src = src.replace(before, after)
    _claim(c, marker)
    _used.append(marker)
    _report.append(f"  patched {marker[:52]:54s} {len(subs)} replacement(s)")
    return nbformat.v4.new_code_cell(src.rstrip())


def md(text):
    return nbformat.v4.new_markdown_cell(text.strip())


def new(text):
    return nbformat.v4.new_code_cell(text.strip())


# ---------------------------------------------------------------- FIXES ----
# (a) submit() declares OUT_FILE_EXACTBYWEIGHT TWICE and never declares the
#     LEDGER #176 above-cap sidecar. The second copy becomes the missing one.
# (b) the wrapper's sidecar delivery loop has no exactByWeightAboveCap entry,
#     so the #176 instrument writes <out>.chrN.exactByWeightAboveCap.txt inside
#     the job and it is DISCARDED. That is why the abovecap reader cell prints
#     "no above-cap sidecar delocalized -- read the count from the job log".
DUP_OUT = ('                "--output", f"OUT_FILE_EXACTBYWEIGHT={od}/{pre}.exactByWeight.txt.gz",\n'
           '                "--output", f"OUT_FILE_EXACTBYWEIGHT={od}/{pre}.exactByWeight.txt.gz",')
NEW_OUT = ('                "--output", f"OUT_FILE_EXACTBYWEIGHT={od}/{pre}.exactByWeight.txt.gz",\n'
           '                "--output", f"OUT_FILE_EXACTBYWEIGHTABOVECAP={od}/{pre}.exactByWeightAboveCap.txt.gz",\n'
           '                "--output", f"OUT_FILE_SPAPINNED={od}/{pre}.spaPinned.txt.gz",')

PAIR_OLD = '"spaFallbacks:${OUT_FILE_SPAFALLBACKS:-}" "exactByWeight:${OUT_FILE_EXACTBYWEIGHT:-}"'
PAIR_NEW = ('"spaFallbacks:${OUT_FILE_SPAFALLBACKS:-}" "exactByWeight:${OUT_FILE_EXACTBYWEIGHT:-}" '
            '"exactByWeightAboveCap:${OUT_FILE_EXACTBYWEIGHTABOVECAP:-}" '
            '"spaPinned:${OUT_FILE_SPAPINNED:-}"')

# (c) the wrapper reads SPA_DUMP but never exports SAIGE_SPA_DUMP_BELOW, which
#     src/SAIGE_test.cpp:131 reads. Without it the LEDGER #172 diagnostic -- dump
#     every CONVERGED saddlepoint below a threshold, the only way to see a
#     confidently-wrong p -- cannot be requested from a dsub job at all.
DUMP_OLD = ('if [[ -n ${SPA_DUMP} ]]; then mkdir -p spadump; '
            'export SAIGE_SPA_DUMP_DIR="${PWD}/spadump"; fi')
DUMP_NEW = (DUMP_OLD + "\n"
            "# LEDGER #172: a CONVERGED saddlepoint below this p is dumped too. An\n"
            "# explicit `if`, which is this file's own stated convention for set -e.\n"
            'if [[ -n ${SPA_DUMP_BELOW:-} ]]; then export SAIGE_SPA_DUMP_BELOW="${SPA_DUMP_BELOW}"; fi')

cells = []

# ================================================================ 0. setup ==
cells += [md("""
# 04b production: SAIGE step 2 on All of Us (FlexRV + SKAT-O), v4 models

Run **sections 0-1 top to bottom on a fresh kernel**; every later cell depends
on them. Section 2 launches production. Section 3 is the LEDGER #163-#176
rerun / splice / check machinery, one parameterised cell per step with
`ANC = "..."` at the top. Section 4 is the reading layer (summaries,
calibration, the discordance screen). Section 5 is LEDGER #175/#176.

**Tidy-up of 2026-10-02.** The one-off investigation cells of 2026-09-29..10-01
are gone; their record is the SAIGE repo (`docs/state/handoff-2026-09-30-rollout-complete.md`,
`handoff-2026-10-01-*.md`, `handoff-2026-10-02-detectability-and-176.md`, LEDGER
#163-#176). Five live defects were fixed in the process:

* **`submit_check_v7` had FOUR definitions.** Two hardcoded the old side to
  `flexrv_v4_fixed/` (trap 14a: that killed 130 of 444 jobs off eas) and one
  fell back to `RUNS[...]["out_prefix"]`, which the 2026-09-30 reruns MUTATED to
  the v7 paths -- so it compared v7 against v7 and passed vacuously. The
  surviving one tries `flexrv_v4_fixed/` then `flexrv_v4/` and SKIPs, naming the
  reason, when neither exists.
* **`classify` meant two things.** Step 4's splice-log verdict and the
  rerun-coverage gap finder. The second is now `rerun_gaps`.
* **`submit()` declared `OUT_FILE_EXACTBYWEIGHT` twice** and never declared the
  LEDGER #176 above-cap sidecar; **the wrapper never delivered it**, so that
  instrument's output was written inside the job and thrown away.
* `stage_check`, `show_log`, `load`, `cct`, `MODEL_DIR`, `submit_summary_v7`
  and `submit_summary_prefix` each existed twice.
* **`submit_abovecap` never passed `geno`**, so it took `submit()`'s `bed`
  default while the coding subsets exist only as `cod.{pgen,pvar,psam}`. That
  is the trap recorded in `54c8d841`: ~61 s, exit 1, reads like localization.
  It now defaults to `pgen`.

**After running this, re-run the `%%writefile` wrapper cell and the upload cell**:
the wrapper changed.

Egress rule: read files into the VM freely, display summaries only, never
print a table of per-cell rows, and write no AoU count below 20 anywhere.
"""),
  md("## 0. Setup: environment, models, traits"),
  old('WORKSPACE_BUCKET = "gs://'),
  old('_S1 = f"{BUCKET}/data/saige/run"'),
  old("TRAITS_BINARY = {"),
  md("### The dsub wrapper and the step-2 wrapper script (verbatim; edit ONLY here, then run the upload cell)"),
  old("%%writefile /home/jupyter/aou_dsub.bash"),
  patched("%%writefile saige_step2_brava_main_original.sh",
          [(PAIR_OLD, PAIR_NEW), (DUMP_OLD, DUMP_NEW)]),
  old("the local wrapper is STALE"),

  md("## 1. Helpers (submit, registry, logs, comparison, paths)"),
  patched("def submit(name, image, masks, fasttest, outdir, anc, pheno,", [(DUP_OUT, NEW_OUT)]),
  old("def build_panel(anc):"),
  md("### Integrity asserts: a failure here stops Run All, which is the point"),
  old("wrapper, helpers and config OK"),
  new('''
# The two 2026-10-02 repairs, asserted so they cannot silently regress.
import inspect
_w = open("saige_step2_brava_main_original.sh").read()
for _needle, _why in (
        ("exactByWeightAboveCap:${OUT_FILE_EXACTBYWEIGHTABOVECAP:-}", "LEDGER #176 above-cap"),
        ("spaPinned:${OUT_FILE_SPAPINNED:-}",                         "LEDGER #174 pinned-samples"),
        ("export SAIGE_SPA_DUMP_BELOW=",                              "LEDGER #172 converged-dump")):
    assert _needle in _w, (
        f"the wrapper does not carry the {_why} line -- it is written inside the job "
        "and discarded, or cannot be requested at all. Re-run the %%writefile cell.")
_s = inspect.getsource(submit)
assert _s.count("OUT_FILE_EXACTBYWEIGHT=") == 1, \\
    "submit() declares OUT_FILE_EXACTBYWEIGHT more than once"
for _needle in ("OUT_FILE_EXACTBYWEIGHTABOVECAP=", "OUT_FILE_SPAPINNED="):
    assert _needle in _s, f"submit() does not declare {_needle}"
print("sidecar delivery OK: #176 above-cap, #174 pinned, #172 converged-dump")
'''),
  md("""### Submission guard: nothing reaches dsub until this cell ARMS it.

Seven cells in this notebook submit jobs when run -- the production launcher, Steps
1, 3, 6, 7, 8 and 10 -- so a Run All on a fresh kernel would launch hundreds. The
guard sits at the one choke point every submitter shares (`aou_dsub` inside a
`subprocess.run`), not on `submit()`, because `submit_splice`, `submit_check_v7`,
`submit_dumpcheck`, `submit_summary*` and `submit_mu_*` all build their own dsub
command and never call `submit()`. **To launch: set `ARMED = True` here, run this
cell, run the step, then set it back.**"""),
  new('''
# False  = nothing submits
# "pre-" = only job names starting with that prefix may submit
# True   = anything may submit (say so on purpose)
#
# A PREFIX is the setting to use. The production launcher cell carries its own
# switch, LAUNCH, baked into its source -- so arming globally for a one-gene
# job and then running that cell is a mass submission with nothing in between.
# ARMED = "fxcap-" lets the one-gene jobs through and still refuses fxv4-/fxv7-.
ARMED = False

import subprocess as _sp
if not getattr(_sp.run, "_dsub_guarded", False):      # idempotent: re-running is safe
    _run_raw = _sp.run
    def _guarded_run(*a, **k):
        cmd = a[0] if a else k.get("args", "")
        probe = " ".join(map(str, cmd)) if isinstance(cmd, (list, tuple)) else str(cmd)
        if "aou_dsub" in probe:
            _m = __import__("re").search(r"--name\\s+(\\S+)", probe)
            _job = _m.group(1).strip("'" + chr(34)) if _m else "<name not found>"
            _ok = (ARMED is True) or (isinstance(ARMED, str) and ARMED and _job.startswith(ARMED))
            if not _ok:
                raise RuntimeError(
                    f"DISARMED: would submit dsub job {_job!r}, but ARMED = {ARMED!r}. "
                    "Set ARMED to True or to a job-name prefix in the submission-guard "
                    "cell, run it, then re-run this one. (subprocess.run was not called.)")
        return _run_raw(*a, **k)
    _guarded_run._dsub_guarded = True
    _sp.run = _guarded_run
    subprocess = _sp            # the notebook calls subprocess.run by that name

def _guard_selftest():
    """Both directions, because a guard that never fires and a guard that always
    fires look the same from a clean notebook."""
    global ARMED
    was = ARMED
    try:
        def _try(label, setting, cmd, want_refusal):
            global ARMED
            ARMED = setting
            try:
                subprocess.run(["bash", "-c", cmd], capture_output=True)
                print(("FAIL " if want_refusal else "ok   ") + label)
            except RuntimeError:
                print(("ok   " if want_refusal else "FAIL ") + label)
        _try("disarmed refuses a submission", False,
             "aou_dsub --name fxcap-x --image i", True)
        _try("prefix arms its own jobs", "fxcap-",
             "aou_dsub --name fxcap-x --image i", False)
        _try("prefix still refuses OTHER jobs", "fxcap-",
             "aou_dsub --name fxv7-amr-x --image i", True)
        ARMED = False
        r = subprocess.run(["bash", "-lc", "echo not-a-submission"], capture_output=True, text=True)
        print("ok   ordinary command still runs:", r.stdout.strip())
    finally:
        ARMED = was
_guard_selftest()
print("ARMED =", ARMED)
'''),
  md("### Constants for the current state (images, scan stamps) and the registry"),
  old('STAMPS = {"sas":'),
  md("### Model QC: variance ratios of every v4 model (flags the 6.91 shape)"),
  old("def vr_ratios_local(anc, tr):"),
  md("### Chromosome bins per ancestry, balanced by group-file size"),
  old("def balance(k):"),
]

# =================================================== 2. production launch ==
cells += [md("""
## 2. Production launch (v4): FlexRV and SKAT-O arms, every model, every ancestry

`RUN_ANCS`, `SPOT` and `LAUNCH` at the top of the launch cell. `already()` keys
on DELIVERED results, so the cell is re-runnable as a mop-up after preemptions;
it refuses to launch while jobs for the same ancestries are RUNNING.
"""),
  # LAUNCH is a SECOND switch, independent of ARMED and baked into this cell's
  # source as True. Start it off: arming for a one-gene job and then running
  # this cell must not be a mass submission.
  patched("RUN_ANCS = (", [("LAUNCH   = True", "LAUNCH   = False")]),
  md("### Monitor: dstat histogram by status, failures by job prefix"),
  old("FAILURE by prefix:"),
  md("### Failure detail (preempted / OOM / other) for the last hours"),
  old("def detail(r):"),
  md("### Delivered-run coverage: chrStatus, partial or refused runs, missing chromosomes"),
  old("trait-arms missing chromosomes"),
]

# ============================================ 3. rerun, splice and check ===
cells += [md("""
## 3. Saddlepoint rerun, splice and check (LEDGER #163-#172)

Steps 1-5 are the v6 rerun-and-splice over the scan's failed regions
(`flexrv_v6/` -> `flexrv_v4_fixed/`), completed for all five ancestries on
2026-09-28 and kept for the next image. Steps 6-7 are the #172 dump. Steps 8-10
are the v7 unfiltered rerun (`flexrv_v7/`) and its control
(`flexrv_v7_checked/`). `ANC = "..."` at the top of each step.
"""),
  md("### Rerun and splice helpers"),
  old("def regions_for(anc, pre, stamp):"),
  old('def submit_splice(v4name, stamp, new_dir='),
  old("def classify(log):"),
  new('''
# The scan stamp is PER ANCESTRY and a wrong one reads as CLEAN (trap 14).
# Derive it from the bucket; never from a scalar.
def stamp_for(anc, arm="flexrv_v4"):
    p = f"{BUCKET}/SAIGE_SLIM_TESTING/spa_scan/{anc}/{arm}/"
    out = subprocess.run(["bash", "-lc", f"gcloud storage ls '{p}'"],
                         capture_output=True, text=True).stdout.split()
    st = sorted(x.rstrip("/").rsplit("/", 1)[-1] for x in out if x.startswith("gs://"))
    if len(st) != 1:
        print(f"  WARNING {anc}: {len(st)} stamps {st} -- using the latest")
    assert st, f"{anc}: no scan stamp under {p} -- any check against it is vacuous"
    return st[-1]
'''),
  md("### Step 1: rerun the binary traits of one ancestry on IMG7"),
  old("todo = [(n, r) for n, r in list(RUNS.items())"),
  md("### Step 2: stage check -- expect `fallbacks 0`, no `CHROMOSOME FAILED`, `STAGE OK`"),
  old('d = f"stage_{ANC}"; os.makedirs'),
  md("### Step 3: splice every rerun job of the ancestry over its v4 files"),
  old("submit_splice(name, STAMPS[ANC])"),
  md("### Step 4: splice verdicts and non-fatal signals"),
  old("kinds[name] = classify(fetch_log("),
  md("### Step 5: positive and negative control over the scan's most extreme v4 burden tails"),
  old("tail_regions.tsv"),
  md("""### Step 6 (LEDGER #172): dump the CONVERGED saddlepoint of the impossible cells.
Needs the wrapper's `SPA_DUMP` block to carry `SAIGE_SPA_DUMP_BELOW` and the repo's
current `aou_spa_dumpcheck.sh` in the home directory. One chromosome per gene; a wrong
one dies at once with `REGIONS matched no region`."""),
  old("BELOW = ["),
  md("""### Step 7: replay each delivered dump inside the SAIGE image, summary lines only.
`carriers` and `raw vector` carry small counts and are left out of the print on purpose."""),
  old("def submit_dumpcheck(anc, tgz, image=IMG8):"),
  old("def replay_lines(name, wait_min=40):"),
  md("""### Step 8 (LEDGER #172 FIXED): rerun every binary FlexRV job on IMG9, unfiltered, into `flexrv_v7/`"""),
  old("def submit_v7(v4name, image=IMG9, spot=False):"),
  old("def submit_rerun_v7(v4name, image=IMG9, spot=False, dry=False):"),
  defonly("def missing_v7(anc):", stop='for a in ["afr", "amr", "eur"]:'),
  defonly("def submitted_names(age=", stop="want = {k.replace("),
  md("""### Step 9: completeness, then the v7 stage check.
`missing_v7` lists jobs whose nine delivered objects are not all present; `stage_check`
runs `aou_spa_stagecheck.sh` over the staged sidecars."""),
  defonly("def stage_check(anc):", stop='for a in ["afr", "amr", "eur"]:', allow_dupes=True),
  new('''
# Per-ancestry accounting from the jobs' own logs. The five counters a step-2
# run prints are the only place #174/#175/#176 surface, and NOTHING in
# harness/tools/ scrapes any of them (checked 2026-10-02) -- so they are
# scraped here until `aou_log_counters.sh` exists.
import re
COUNTERS = {
    "fallbacks":  r"SPA fallbacks: (\\d+) ",
    "exact_wt":   r"Exact test by weight: (\\d+) cell",
    "above_cap":  r"Concentrated cells above the exact test's enumeration cap: (\\d+) cell",
    "pinned":     r"Pinned samples: (\\d+) ",
}
def counters(anc, prefix="fxv7-"):
    names = sorted(n for n in RUNS if n.startswith(f"{prefix}{anc}-"))
    tot = {k: 0 for k in COUNTERS}
    for n in names:
        log = fetch_log(n)
        per = {k: sum(int(x) for x in re.findall(p, log)) for k, p in COUNTERS.items()}
        neff = re.findall(r"Effective sample size \\([^)]*\\): ([\\d.]+) of (\\d+)", log)
        frac = f"{100*float(neff[0][0])/float(neff[0][1]):.1f}%" if neff else "-"
        print(f"  {n:38s} " + "  ".join(f"{k} {per[k]:>5d}" for k in COUNTERS)
              + f"   Neff {frac}")
        for k in COUNTERS: tot[k] += per[k]
    print(f"{anc}: {len(names)} jobs | " + " | ".join(f"{k} {v}" for k, v in tot.items()))
    return tot

ANC = "sas"
stage_check(ANC)
counters(ANC)
'''),
  md("""### Step 10: the control for a COMPLETE rerun (`REGIONS=ALL FAILURES=NONE`).
The old side is the current final -- the spliced copy if there is one, else v4. It is
resolved by LOOKING, not assumed: a hardcoded `flexrv_v4_fixed/` passes the whole eas
chain and killed 130 of 444 jobs elsewhere (trap 14a), and a fallback to
`RUNS[...]["out_prefix"]` compares v7 against v7, because the 2026-09-30 reruns mutated
that field to the v7 paths."""),
  defonly("SKIP -- no old side in flexrv_v4_fixed or flexrv_v4", stop="import collections"),
  new('''
ANC = "sas"
for k in sorted(n for n in RUNS if n.startswith(f"fxv4-{ANC}-") and RUNS[n]["trait"] == "binary"):
    p7 = RUNS[k]["out_prefix"].rsplit("/", 1)[-1].replace("fxv4", "fxv7")
    submit_check_v7(p7, anc=ANC)
'''),
  new('''
# The control's verdict lines, no per-cell rows.
ANC = "sas"
for name in sorted(n.replace("fxv4", "check7") for n in RUNS
                   if n.startswith(f"fxv4-{ANC}-") and RUNS[n]["trait"] == "binary"):
    log = fetch_log(name).splitlines()
    keep = [l for l in log if any(k in l for k in
            ("ERROR", "WARNING", "REGIONS=ALL", "cells routed to the exact test",
             "unobservable burden cells", "changed row(s)", "spliced", "carried over",
             "pseudo-marker rows explained"))]
    print("=====", name, "| ROWS lines:", sum(l.startswith("ROWS ") for l in log))
    print("\\n".join(l[:170] for l in keep))
'''),
  md("### The one-off check of named prefixes (the four blocked ones of 2026-10-01)"),
  old("def submit_check_one(anc, pre, image="),
  md("### Rerun-coverage gaps: jobs with no fixed copy, split into benign and real"),
  patched("def classify(anc):", [("def classify(anc):", "def rerun_gaps(anc):"),
                                 ("    classify(a)", "    rerun_gaps(a)")]),
]

# =============================================== 4. summaries and reading ==
cells += [md("""
## 4. Summaries, calibration and the FlexRV-vs-SKAT-O screen

`summary_v4/` is built from `flexrv_v7_checked/` when present; `summary_v4_prefix/`
forces the pre-fix fallback with a non-existent `FLEXRV_V7_DIR`, which is how the
old-vs-new comparison gets both sides. **The log's `K v7` count is the only thing
that distinguishes a current table from a stale one -- the filename does not**
(trap 14c), so the provenance cell below is not optional.
"""),
  md("### Submitters: post-fix, pre-fix, and the v7-directory variant"),
  defonly("def submit_summary(anc, image=", stop='submit_summary("eas")'),
  defonly("def submit_summary_v7(anc, image=", stop='submit_summary_v7("sas")',
          allow_dupes=True),
  defonly("__does_not_exist__", stop='for anc in ["sas", "amr", "eur", "afr"]:'),
  md("### Provenance: what each summary log says it was built from (trap 14c)"),
  new('''
# The `flexrv: N v4 files, M fixed, K v7` line is the provenance. K is derived
# from RUNS, not pinned: the hardcoded (sas 6, amr 122, eur 320) went stale.
def summary_provenance(ancs=("sas", "amr", "afr", "eur", "eas"), job="summary"):
    for anc in ancs:
        want = len({RUNS[k]["out_prefix"].rsplit("/", 1)[-1].replace("fxv4", "fxv7")
                    for k in RUNS if k.startswith(f"fxv4-{anc}-")
                    and RUNS[k]["trait"] == "binary"})
        p = log_path(f"{job}-{anc}")
        if not p:
            print(f"{anc}: no {job} log yet (expect {want} v7)"); continue
        date = p.split("/")[-3] if len(p.split("/")) > 3 else "?"
        txt = subprocess.run(["bash", "-lc", f"gcloud storage cat '{p}'"],
                             capture_output=True, text=True).stdout
        line = next((l for l in txt.splitlines() if "flexrv:" in l), "(no flexrv line)")
        got = re.search(r"(\\d+) v7", line)
        ok = "OK " if got and int(got.group(1)) == want else "STALE?"
        print(f"{ok} {anc} (expect {want} v7)  log dated {date}\\n     {line.split('] ')[-1]}")
summary_provenance()
'''),
  md("### Per-ancestry pipeline coverage: reruns, splices, summaries"),
  new('''
def coverage(ancs=("afr", "amr", "eur", "sas", "eas")):
    def n(path, pat):
        return int(subprocess.run(["bash", "-lc",
            f"gcloud storage ls '{BUCKET}/SAIGE_SLIM_TESTING/{path}/{pat}' 2>/dev/null | wc -l"],
            capture_output=True, text=True).stdout.strip() or 0)
    for anc in ancs:
        want = len([k for k in RUNS if k.startswith(f"fxv4-{anc}-")
                    and RUNS[k]["trait"] == "binary"])
        print(f"{anc}: reruns {n(f'{anc}/flexrv_v7', '*.chrStatus.txt.gz')}/{want} | "
              f"splices {n(f'{anc}/flexrv_v7_checked', '*.splice.txt')}/{want} | "
              f"summary {'yes' if n('summary_v4', f'{anc}.flexrv_vs_skato.v4.tsv.gz') else 'no'} | "
              f"old-vs-new {'yes' if n('summary_v4_prefix', f'{anc}.flexrv_vs_skato.v4.tsv.gz') else 'no'}")
coverage()
'''),
  md("### Log reader (filtered) and the delivered-object audit"),
  defonly("def log_path(job_glob):", stop='show_log("check7-'),
  new('''
# Which sidecars were actually delivered, and whether exactByWeight is a real
# file or the wrapper's placeholder.
import collections
SFX = ("exactByWeightAboveCap.txt.gz", "exactByWeight.txt.gz", "spaFallbacks.txt.gz",
       "txt.gz", "singleAssoc.txt.gz", "pooledTests.txt.gz")
def delivered_audit(anc, arm="flexrv_v7"):
    out = subprocess.run(["bash", "-lc",
        f"gcloud storage ls -l '{BUCKET}/SAIGE_SLIM_TESTING/{anc}/{arm}/*'"],
        capture_output=True, text=True).stdout.splitlines()
    size, have = {}, collections.defaultdict(set)
    for line in out:
        f = line.split()
        if len(f) < 3 or not f[-1].startswith("gs://"): continue
        name = f[-1].rsplit("/", 1)[-1]
        for s in SFX:                       # longest suffixes first, so they win
            if name.endswith("." + s):
                pre = name[: -(len(s) + 1)]; have[pre].add(s); size[(pre, s)] = int(f[0]); break
    for pre in sorted(have):
        ew = "exactByWeight.txt.gz"; n = size.get((pre, ew))
        v = "MISSING" if ew not in have[pre] else ("PLACEHOLDER" if n < 200 else f"{n} bytes")
        print(f"{pre:<44} exactByWeight: {v:<16} objects: {len(have[pre])}")
delivered_audit("eas")
'''),
  md("### Old vs new: did the #172 fix move what it should and nothing else?"),
  new('''
# One parameterised control replacing the five per-ancestry copies. The assert
# is the point: without it a pair of tables that are both post-fix compares
# clean and says nothing.
import pandas as pd, numpy as np
def oldnew(anc):
    def get(d):
        f = f"/tmp/{anc}.{d}.tsv.gz"
        subprocess.run(["bash", "-lc", f"gcloud storage cp "
            f"'{BUCKET}/SAIGE_SLIM_TESTING/{d}/{anc}.flexrv_vs_skato.v4.tsv.gz' {f}"], check=True)
        return pd.read_csv(f, sep="\\t")
    old_t, new_t = get("summary_v4_prefix"), get("summary_v4")
    assert "v7" in set(new_t.source) and "v7" not in set(old_t.source), (
        f"{anc}: these are not pre-fix vs post-fix "
        f"(old {sorted(set(old_t.source))}, new {sorted(set(new_t.source))})")
    j = old_t.merge(new_t, on=["anc", "trait", "gene"], suffixes=("_old", "_new"))
    j = j[j.p_flexrv_old.notna() & j.p_flexrv_new.notna()].copy()
    b = j[j.source_new == "v7"].copy(); b["ratio"] = b.p_flexrv_new / b.p_flexrv_old
    q = j[j.source_new != "v7"]
    print(f"{anc}: binary genes {len(b)} | unchanged {(b.ratio == 1).mean():.4f} | "
          f">2x up {int((b.ratio > 2).sum())} | <0.5x {int((b.ratio < 0.5).sum())} | "
          f"1e-8 -> >1e-3 {int(((b.p_flexrv_old < 1e-8) & (b.p_flexrv_new > 1e-3)).sum())}")
    print(f"     quantitative (not rerun): {len(q)} genes, "
          f"{int((q.p_flexrv_old != q.p_flexrv_new).sum())} changed  <- must be 0")
    # WHICH genes moved, not just how many: the counts say the fix behaved, the
    # names are what you actually look at afterwards.
    top = b.nlargest(10, "ratio")
    print(top[["trait", "gene", "p_flexrv_old", "p_flexrv_new",
               "set_best_new", "ratio"]].to_string(index=False))
    return b
for a in ["sas", "amr", "afr", "eur"]:
    oldnew(a)
'''),
  md("### Per-trait null calibration from the summary table"),
  new('''
import math
def calibration(anc, source="v7"):
    d = pd.read_csv(f"/tmp/{anc}.summary_v4.tsv.gz", sep="\\t")
    b = d[(d.source == source) & d.p_flexrv.notna()]
    rows = []
    for tr, t in b.groupby("trait"):
        n = len(t); obs = t.p_flexrv.min(); med = 1 - 2 ** (-1 / n)
        top = t.nsmallest(1, "p_flexrv").iloc[0]
        rows.append(dict(trait=tr, genes=n, min_p=obs, null_med=med, ratio=obs / med,
                         P_ge=(1 - obs) ** n, set_best=top.set_best, p_skato=top.p_skato,
                         orders=(math.log10(top.p_skato / obs)
                                 if top.p_skato and obs else float("nan"))))
    print(f"=== {anc} ===")
    print(pd.DataFrame(rows).sort_values("ratio")
          .to_string(index=False, float_format=lambda v: f"{v:.3g}"))
calibration("sas")
'''),
  md("""### The FlexRV-vs-SKAT-O discordance screen.
`p_skato` SATURATES -- SKAT-O returns exactly 1 once its smallest per-rho p exceeds
~0.78 (LEDGER #86) -- so the most NULL cells sort to the top as the most suspicious
unless rows above ~0.9 are excluded, which `screen()` does. `ur_only` flags the
MAC1/3/5/10 rungs, which sit at or below the collapse cutoff and therefore test the
pseudo-marker alone: the #163/#164 signature."""),
  new('''
# P is the LOCAL root holding summary_v4/ and summary_v4_prefix/. Set it and
# stage them once; the screen reads nothing from the bucket.
P = os.path.expanduser("~/summ")
os.makedirs(P, exist_ok=True)
for d in ("summary_v4", "summary_v4_prefix"):
    os.makedirs(f"{P}/{d}", exist_ok=True)
    subprocess.run(["bash", "-lc",
        f"gcloud storage cp '{BUCKET}/SAIGE_SLIM_TESTING/{d}/*.tsv.gz' '{P}/{d}/'"],
        capture_output=True)
print(subprocess.run(["bash", "-lc", f"ls -1 {P}/summary_v4 | head"],
                     capture_output=True, text=True).stdout)
'''),
  defonly("EPS = 1e-320", stop='q = screen("eur")'),
  new('''
q = screen("eur")
q = q[q.source != "v7"]                      # quantitative = the artefact-immune arm
top = q.nlargest(20, "o4")
print(top[["trait","gene","p_flexrv","p_tight4","o4","omin","set_best","n_sets","ur_only"]]
        .to_string(index=False))
print("\\nof the top 20:  omin > 10:", int((top.omin > 10).sum()),
      "| pseudo-marker-only set:", int(top.ur_only.sum()),
      "| distinct traits:", top.trait.nunique())
'''),
  md("### Pooled-row provenance: which CCT rule fired (LEDGER #86/#127)"),
  new('''
def pooled_rules(anc, pre, arm="flexrv_v7_checked"):
    base = f"{BUCKET}/SAIGE_SLIM_TESTING/{anc}/{arm}/{pre}"
    txt = subprocess.run(["bash", "-lc",
        f"gcloud storage cat '{base}.pooledTests.txt.gz' | gzip -dc"],
        capture_output=True, text=True).stdout.splitlines()
    if not txt:
        print("no pooledTests for", pre); return
    hdr = txt[0].split("\\t"); rows = [l.split("\\t") for l in txt[1:] if l.strip()]
    print(f"{len(rows)} disclosed regions | columns: {hdr}")
    ri = next((i for i, h in enumerate(hdr) if "rule" in h.lower() or "cct" in h.lower()), -1)
    if ri >= 0:
        for rule, n in collections.Counter(r[ri] for r in rows).most_common():
            print(f"  {n:6d}  {rule}")
'''),
]

# ================================================ 5. the mu / #175 / #176 ==
cells += [md("""
## 5. Fitted-value decades, pinned samples, and LEDGER #175 / #176

`aou_mu_decades.sh` and `aou_mu_drivers.sh` characterise the fitted-value
distribution of every model; `aou_pinned_scan.sh` counts pinned samples
(LEDGER #174), masked below 20. `submit_abovecap` runs one gene on IMG10,
which reports the #176 above-cap count and #175's effective sample size.
"""),
  defonly("def submit_mu_decades(anc, image=IMG9):", stop='for anc in ["sas", "afr", "eur", "amr"]:'),
  defonly("def submit_mu_drivers(anc, traits, image=IMG9):", stop='submit_mu_drivers("eur",'),
  defonly("the notebook's STEP1:", stop='for anc in ["eas", "sas", "afr", "eur", "amr"]:'),
  md("### Read the decade tables and the pinned scan"),
  new('''
DEC = ["le1e-15","1e-15_1e-12","1e-12_1e-10","1e-10_1e-9","1e-9_1e-8",
       "1e-8_1e-7","1e-7_1e-6","1e-6_1e-5"]
def mu_decades(ancs=("eas","sas","afr","eur","amr")):
    tabs = []
    for anc in ancs:
        f = f"/tmp/mu_decades.{anc}.tsv.gz"
        if subprocess.run(["bash", "-lc", f"gcloud storage cp "
            f"'{BUCKET}/SAIGE_SLIM_TESTING/mu_decades/mu_decades.{anc}.tsv.gz' {f}"],
            capture_output=True).returncode == 0:
            t = pd.read_csv(f, sep="\\t", dtype=str); t["anc"] = anc; tabs.append(t)
    d = pd.concat(tabs, ignore_index=True)
    d["trait"] = d.model.str.replace(r"^[a-z]+_|\\.rda$", "", regex=True)
    flag = d[d["n_le_1e-06"] != "0"]
    print(f"{len(d)} binary models, {len(flag)} flagged at 1e-6\\n")
    print(flag[["anc","trait","N","cases"] + DEC]
          .sort_values(["trait","anc"]).to_string(index=False))
    return d, flag

def pinned_scan():
    t = load(f"{BUCKET}/SAIGE_SLIM_TESTING/pinned/pinned.v4.tsv.gz")
    print("models:", len(t))
    print(t.status.value_counts().to_string()); print(t.gates.value_counts().to_string())
    # the two controls FIRST: eas/MatHem is the model known to carry pinned
    # samples, sas is the arm that should carry none. Read them before the
    # table -- a scan that reports nothing anywhere is indistinguishable from
    # a scan that did not run.
    print("\\npositive control (eas/MatHem):")
    print(t[(t.anc == "eas") & (t.trait == "mathem")].to_string(index=False))
    print("\\nnegative control (sas status counts):")
    print(t[t.anc == "sas"].status.value_counts().to_string())
    w = t[t.status == "warning"].sort_values(["anc", "trait"])
    print(w[["anc","trait","pinned","N","percent","min_fitted","max_fitted"]].to_string(index=False))
    return t
'''),
  md("""### LEDGER #176: one gene on IMG10, which counts the concentrated cells the exact
test could not enumerate. **Observation only** -- no p-value and no routing decision
moves. With the wrapper repair above, the `.exactByWeightAboveCap.txt.gz` sidecar is
now DELIVERED instead of being written inside the job and discarded."""),
  # LEDGER/recipe trap (54c8d841, 2026-10-02): submit() defaults geno to "bed",
  # but the coding subsets exist only as cod.{pgen,pvar,psam}. A caller that
  # omits geno dies in ~61 s with exit 1, which reads like a localization
  # problem. submit_abovecap omitted it; it now defaults to pgen itself.
  patched("def submit_abovecap(anc, pheno, gene, image=IMG10",
          [('trait="binary", **kw):', 'trait="binary", geno="pgen", extra_env=None, **kw):'),
           ("                  pooled=True, sidecars=True,",
            "                  geno=geno, pooled=True, sidecars=True,"),
           ('"RCORR": "1", "REGIONS": gene}, **kw)',
            '"RCORR": "1", "REGIONS": gene, **(extra_env or {})}, **kw)')]),
  new('''
# Read one above-cap run: the sidecar's concentration numbers, then the cells
# themselves. n_material is the number that decides whether raising
# kExactByWeightCarriers (13 today, 2**13 enumerations) could ever reach them.
import io, gzip, csv
def read_abovecap(anc, pre, outdir=None):
    OUT = outdir or f"{BUCKET}/SAIGE_SLIM_TESTING/{anc}/flexrv_abovecap"
    def cat(p):
        r = subprocess.run(["gcloud","storage","cat",p], capture_output=True)
        if r.returncode != 0: return None
        return gzip.decompress(r.stdout).decode() if r.stdout[:2] == b"\\x1f\\x8b" else r.stdout.decode()
    side = cat(f"{OUT}/{pre}.exactByWeightAboveCap.txt.gz")
    if not side:
        print("!! no above-cap sidecar delocalized. If this job predates the wrapper "
              "repair of 2026-10-02 the file was written inside the job and discarded; "
              "read the count from the log instead."); return None
    rows = list(csv.DictReader(io.StringIO(side), delimiter="\\t"))
    print(f"=== {pre}: {len(rows)} above-cap cells "
          "(n_eff <= --max_MAC_for_ER, n_material > 13) ===")
    for r in rows:
        print(f"  {r['Weight']:26s} n_eff {float(r['n_eff']):7.3f}  "
              f"n_material {r['n_material']:>4s}  MAC {r['MAC']}")
    mat = sorted(int(r["n_material"]) for r in rows)
    if mat:
        print(f"n_material: min {mat[0]}, median {mat[len(mat)//2]}, max {mat[-1]} "
              f"(cap 13 = 2**13 enumerations; 2**20 ~ 1e6)")
    return rows
'''),
  new("""
# The whole #176 replay in one call. Four things it does that were four manual
# steps, each of which can fail quietly:
#   - builds the tarball path from the prefix, so it cannot be mistyped;
#   - REFUSES an empty tar. dsub always delivers the declared output, so a run
#     where nothing fell below SPA_DUMP_BELOW looks identical to a successful
#     one until the replay reports `dumps: 0`;
#   - defaults to the image that PRODUCED the dumps, not IMG8. aou_spa_dumpcheck
#     replays the run's OWN kernel, so an older image compares the vector
#     against different arithmetic and can come back clean for the wrong reason;
#   - names the job from the prefix rather than split('-')[3], which was written
#     for the fxbelow-<anc>-<trait>-<gene> shape and yields '22' off chr1-22.
def dumpcheck_abovecap(anc, pre, image=IMG10, outdir="flexrv_abovecap", wait_min=40):
    tgz = f"{BUCKET}/SAIGE_SLIM_TESTING/{anc}/{outdir}/{pre}.spaDump.tgz"
    out = subprocess.run(["bash", "-lc", f"gcloud storage ls -l '{tgz}'"],
                         capture_output=True, text=True).stdout.split()
    size = int(out[0]) if out and out[0].isdigit() else 0
    if size < 1000:
        print(f"NO DUMPS: {tgz.rsplit('/', 1)[-1]} is {size} bytes. Nothing fell below "
              "SPA_DUMP_BELOW on that run -- raise the threshold and rerun, do not "
              "read this as 'the saddlepoint is fine'.")
        return None
    name = f"dumpcheck-{anc}-{pre.split('.')[1].lower()}"[:63]
    args = ["aou_dsub", "--name", name, "--image", image, "--input", f"IN_DUMP={tgz}",
            "--machine-type", "n2-standard-2", "--disk-size", "20",
            "--script", "aou_spa_dumpcheck.sh"]
    script = "source /home/jupyter/aou_dsub.bash\\n" + " ".join(shlex.quote(a) for a in args)
    r = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=os.environ)
    print(f"{name}  ({size} bytes)  |", (r.stdout or r.stderr).strip().splitlines()[-1])
    replay_lines(name, wait_min)
    return name

# ARMED = True; dumpcheck_abovecap("amr", "fxcap-amr-efrmb.EFRMB.amr.chr1-22"); ARMED = False
"""),
  md("### Utilities: a compact dstat. `tail_log`, `show_log` and `fetch_log` are\nalready defined in section 1."),
  new('''
get_ipython().system('dstat --provider google-batch --project "${GOOGLE_CLOUD_PROJECT}" '
                     '--users \\'*\\' --status \\'*\\' --age 2h')
'''),
]

_bad = 0
for _i, _c in enumerate(cells):
    if _c.cell_type != "code" or (_c.source or "").lstrip().startswith(("%%", "!")):
        continue
    try:
        compile(_c.source, f"cell{_i}", "exec")
    except SyntaxError as _e:
        _bad += 1
        print(f"SYNTAX ERROR in generated cell {_i} line {_e.lineno}: {_e.msg}")
        print(f"   {(_e.text or '').rstrip()[:100]}")
if _bad:
    sys.exit(f"{_bad} generated cell(s) do not compile -- nothing written")

out = nbformat.v4.new_notebook(cells=cells, metadata=nb.metadata)
nbformat.write(out, NEW)
ncode = sum(1 for c in cells if c.cell_type == "code")
print("\n".join(_report))
print(f"\nwrote {NEW}: {len(cells)} cells ({ncode} code, {len(cells)-ncode} markdown) "
      f"from {len(nb.cells)} cells ({len(code)} code) in {OLD}")
# The builder's OWN drop list: exact, because it is the complement of what was
# claimed. A check done afterwards by comparing sources cannot see a dropped
# cell that is byte-identical to a kept one.
dropped = [(i, c) for i, c in enumerate(nb.cells)
           if c.cell_type == "code" and id(c) not in _seen]
print(f"kept {len(_used)} code cells, dropped {len(dropped)}:")
for i, c in dropped:
    head = next((l.strip() for l in (c.source or "").splitlines()
                 if l.strip() and not l.strip().startswith("#")), "(empty)")
    twin = any(id(k) in _seen and k.source.rstrip() == c.source.rstrip()
               for k in nb.cells if k.cell_type == "code")
    print(f"  {i:>4}{'  [identical twin KEPT]' if twin else '                      '}  {head[:88]}")
