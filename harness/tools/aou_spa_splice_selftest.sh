#!/bin/bash
# Self-test for aou_spa_splice.sh and aou_spa_stagecheck.sh on FABRICATED
# files (no AoU data): one positive control that must splice, and planted
# defects that must each be refused for the reason planted. Run from anywhere:
#   bash harness/tools/aou_spa_splice_selftest.sh
# The fixtures encode every structure the 2026-09-27/28 rollout met:
# pseudo-marker and real-marker singleAssoc rows, a marker shared by
# overlapping regions (multiset), a burden cell with Number_rare > 0
# (unobservable in v4), pseudo-marker rows changed by their own exact-test
# eligibility (LEDGER #166) bounded by the 0/1-weight twin, a pooled sidecar
# placeholder, and a chrStatus with a failed chromosome.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d); trap 'rm -rf "${T}"' EXIT
cd "${T}"
python3 - <<'PY'
import gzip, os
def w(path, lines, gz=True):
    data = "\n".join(lines) + "\n"
    (gzip.open(path, "wt") if gz else open(path, "w")).write(data)
H = "Region\tGroup\tmax_MAF\tWeight\tPvalue_Burden\tBETA_Burden\tSE_Burden\tMAC\tMAC_case\tMAC_control\tNumber_rare\tNumber_ultra_rare"
def row(reg, w_, p, mac, nr=0): return f"{reg}\tgrp\t0.001\t{w_}\t{p}\t1\t1\t{mac}\t1\t{mac-1}\t{nr}\t{mac-nr}"
def cau(reg, p): return f"{reg}\tCauchy\tNA\tNA\t{p}\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
SH = "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\tMissingRate\tBETA\tSE\tTstat\tvar\tp.value\tp.value.NA\tIs.SPA"
def ur(reg, w_, ac, p, na, spa): return f"UR\tUR\t{reg}:grp:0.001:{w_}\tUR\tUR\t{ac}\t0.0001\t0\t1\t1\t1\t1\t{p}\t{na}\t{spa}"
def rm(pos, ac, p, na, spa): return f"1\t{pos}\tchr1:{pos}:A:G\tA\tG\t{ac}\t0.0005\t0\t1\t1\t1\t1\t{p}\t{na}\t{spa}"
# regions: R1 (recorded failure W1, MAC-3 cell W2, unchanged W3, unobservable W4), R2 (nothing changes),
# R4 (tail region: only an unobservable cell changes), R5 (unw_unw count 3: weighted set's pseudo-marker changes),
# R6 (lof_unw count 2: lof x mafbeta changes; cbrt x unw must NOT change), R7 (lof_unw count 9), R3 (not rerun)
old_res = [H, row("R1","W1","1e-100",9), row("R1","W2","0.02",3), row("R1","W3","0.4",9), row("R1","W4","1.6e-7",47,2), cau("R1","1e-99"),
           row("R2","W1","0.3",7), cau("R2","0.3"), row("R4","W1","2e-8",30,5), row("R4","W2","0.5",30,0), cau("R4","4e-8"),
           row("R5","AM_unw_unw","0.2",9,2), row("R5","AM_cbrt0.2_mafbeta2","0.3",9,2), cau("R5","0.25"),
           row("R6","AM_unw_unw","0.2",30,6), row("R6","AM_lof_MAC100","0.3",12,4), row("R6","AM_cbrt0.2_unw","0.4",30,6), row("R6","AM_lof_unw","0.3",12,4), row("R6","AM_lof_mafbeta2","0.3",12,4), cau("R6","0.3"),
           row("R7","AM_unw_unw","0.2",40,10), row("R7","AM_lof_unw","0.3",30,8), row("R7","AM_lof_mafbeta2","0.3",30,8),
           row("R3","W1","0.5",6), cau("R3","0.5")]
new_res = [H, row("R1","W1","0.001",9), row("R1","W2","0.05",3), row("R1","W3","0.4",9), row("R1","W4","0.068",47,2), cau("R1","0.002"),
           row("R2","W1","0.3",7), cau("R2","0.3"), row("R4","W1","0.03",30,5), row("R4","W2","0.5",30,0), cau("R4","0.06"),
           row("R5","AM_unw_unw","0.2",9,2), row("R5","AM_cbrt0.2_mafbeta2","0.3",9,2), cau("R5","0.25"),
           row("R6","AM_unw_unw","0.2",30,6), row("R6","AM_lof_MAC100","0.3",12,4), row("R6","AM_cbrt0.2_unw","0.4",30,6), row("R6","AM_lof_unw","0.3",12,4), row("R6","AM_lof_mafbeta2","0.3",12,4), cau("R6","0.3"),
           row("R7","AM_unw_unw","0.2",40,10), row("R7","AM_lof_unw","0.3",30,8), row("R7","AM_lof_mafbeta2","0.3",30,8)]
shared = rm(500,20,"0.7","0.7","false")   # a marker in three regions: three identical old copies, one rerun copy
old_sa = [SH, ur("R1","W1",5,"1e-100","1e-100","false"), ur("R1","W2",5,"0.02","0.01","true"), ur("R1","W3",5,"0.4","0.4","false"),
          ur("R2","W1",5,"0.3","0.3","false"), ur("R3","W1",5,"0.5","0.5","false"),
          rm(100,12,"1e-90","1e-90","false"), rm(200,3,"0.02","0.01","true"), rm(300,15,"0.3","0.3","false"), rm(400,11,"0.2","0.2","false"),
          shared, shared, shared,
          ur("R5","AM_unw_unw",3,"0.5","0.5","false"), ur("R5","AM_cbrt0.2_mafbeta2","0.61","0.4","0.02","true"),
          ur("R6","AM_unw_unw",9,"0.5","0.5","false"), ur("R6","AM_lof_MAC100",2,"0.3","0.02","true"), ur("R6","AM_cbrt0.2_unw","0.6","0.2","0.02","true"),
          ur("R6","AM_lof_unw",2,"0.3","0.02","true"), ur("R6","AM_lof_mafbeta2","1.98","0.3","0.02","true"),
          ur("R7","AM_unw_unw",30,"0.5","0.5","false"), ur("R7","AM_lof_unw",9,"0.3","0.02","true"), ur("R7","AM_lof_mafbeta2","8.1","0.3","0.02","true")]
new_sa = [SH, ur("R1","W1",5,"0.001","1e-100","true"), ur("R1","W2",5,"0.05","0.01","true"), ur("R1","W3",5,"0.4","0.4","false"),
          ur("R2","W1",5,"0.3","0.3","false"),
          rm(100,12,"0.002","1e-90","true"), rm(200,3,"0.03","0.01","true"), rm(300,15,"0.3","0.3","false"), shared,
          ur("R5","AM_unw_unw",3,"0.5","0.5","false"), ur("R5","AM_cbrt0.2_mafbeta2","0.61","0.45","0.02","true"),
          ur("R6","AM_unw_unw",9,"0.5","0.5","false"), ur("R6","AM_lof_MAC100",2,"0.35","0.02","true"), ur("R6","AM_cbrt0.2_unw","0.6","0.2","0.02","true"),
          ur("R6","AM_lof_unw",2,"0.35","0.02","true"), ur("R6","AM_lof_mafbeta2","1.98","0.36","0.02","true"),
          ur("R7","AM_unw_unw",30,"0.5","0.5","false"), ur("R7","AM_lof_unw",9,"0.3","0.02","true"), ur("R7","AM_lof_mafbeta2","8.1","0.3","0.02","true")]
def variant(sa, drop, add): return [l for l in sa if l != drop] + [add]
plants = {
  "ok":   (new_res, new_sa),
  "real": (new_res, variant(new_sa, rm(300,15,"0.3","0.3","false"), rm(300,15,"0.31","0.3","false"))),          # unexplained real-marker change
  "out":  (new_res, new_sa + [ur("R9","W1",5,"0.1","0.1","false")]),                                               # pseudo-marker row outside the list
  "obs":  ([l if not l.startswith("R2\tgrp\t0.001\tW1\t") else row("R2","W1","0.31",7) for l in new_res], new_sa),  # observable cell changed unexplained
  "dup":  (new_res, variant(new_sa, shared, rm(500,20,"0.71","0.7","false"))),                                     # shared marker's rerun row differs
  "ur":   (new_res, variant(new_sa, ur("R6","AM_cbrt0.2_unw","0.6","0.2","0.02","true"), ur("R6","AM_cbrt0.2_unw","0.6","0.25","0.02","true"))),  # weighted set, maximal count 9
  "twin": (new_res, variant(new_sa, ur("R7","AM_lof_mafbeta2","8.1","0.3","0.02","true"), ur("R7","AM_lof_mafbeta2","8.1","0.31","0.02","true"))),  # lof x mafbeta, twin count 9
  "ew":   ([l if not l.startswith("R2\tgrp\t0.001\tW1\t") else row("R2","W1","0.31",7) for l in new_res], new_sa),  # the "obs" change, listed by the rerun's exactByWeight sidecar (LEDGER #172)
}
# LEDGER #174: the same two changes, explained (or not) by the rerun's spaPinned sidecar
plants["pn"]     = (plants["obs"][0], new_sa)      # the "obs" cell, listed as a guarded burden test: accepted
plants["pnref"]  = (plants["obs"][0], new_sa)      # listed, but the guard REFUSED (over the budget): nothing was left out, so nothing is explained
# LEDGER #176: the same change, explained (or not) by the rerun's stretchGate sidecar
plants["sg"]     = (plants["obs"][0], new_sa)      # the "obs" cell, listed with reported = exact: accepted
plants["sgspa"]  = (plants["obs"][0], new_sa)      # listed (above the floor) but reported = saddlepoint: nothing is explained
plants["sghdr"]  = (plants["obs"][0], new_sa)      # the concatenated sidecar: a 9-column header from a chromosome with no gated cell, the gated row below it in 11 fields
plants["pnreal"] = (new_res, plants["real"][1])    # the "real" marker, listed as a guarded single test: accepted
plants["pnur"]   = (new_res, plants["ur"][1])      # the "ur" pseudo-marker row, listed as a guarded ultra_rare test: accepted
# a COMPLETE rerun (REGIONS=ALL, FAILURES=NONE; LEDGER #172): every region of the old table, one routed cell changed
full_res = [l if not l.startswith("R2\tgrp\t0.001\tW1\t") else row("R2","W1","0.31",7) for l in old_res]
full_res = [l if l != cau("R2","0.3") else cau("R2","0.31") for l in full_res]
full_sa  = [l if l != ur("R2","W1",5,"0.3","0.3","false") else ur("R2","W1",5,"0.31","0.3","true") for l in old_sa]
# CONTROL 1c (2026-09-30): statistic identity. "ur" and "twin" above move p.value
# ALONE, which is the measured signature of a p-value-method change (#166/#172) --
# they are now ACCEPTED by that tier and are its positive controls. These three are
# the negative controls: the same rows with a column the exact test never touches
# moved as well, which must still be refused. `var` is the discriminating one (a
# perturbed variance ratio moves var, BETA and p.value.NA but NOT Tstat).
def urvar(p_, var_="1", beta="1"):  # R6 AM_cbrt0.2_unw with var/BETA overridden
    return f"UR\tUR\tR6:grp:0.001:AM_cbrt0.2_unw\tUR\tUR\t0.6\t0.0001\t0\t{beta}\t1\t1\t{var_}\t{p_}\t0.02\ttrue"
base_ur = ur("R6","AM_cbrt0.2_unw","0.6","0.2","0.02","true")
plants["sidvar"]  = (new_res, variant(new_sa, base_ur, urvar("0.25", var_="1.001")))   # var moved too
plants["sidbeta"] = (new_res, variant(new_sa, base_ur, urvar("0.25", beta="1.001")))   # BETA moved too
plants["sidac"]   = (new_res, variant(new_sa, base_ur,                                  # AC_Allele2 moved too
    "UR\tUR\tR6:grp:0.001:AM_cbrt0.2_unw\tUR\tUR\t0.61\t0.0001\t0\t1\t1\t1\t1\t0.25\t0.02\ttrue"))
plants["allok"]  = (full_res, full_sa)   # the sidecar lists R2|W1: accepted
plants["allbad"] = (full_res, full_sa)   # no sidecar: the same change is unexplained
# 2026-10-08: a complete rerun made with SINGLE_IN_GROUP=FALSE delivers the wrapper's one-line placeholder
# instead of a single-variant file: the results splice as allok, the placeholder is NOTED and no
# single-variant file is written (the old one's pseudo-marker rows would be stale for every rerun region)
plants["nosingle"] = (full_res, ["# no singleAssoc output: SINGLE_IN_GROUP=FALSE (variant results come from the SKAT-O run)"])
os.makedirs("in", exist_ok=True)
w("in/fxv4-test.txt.gz", old_res); w("in/fxv4-test.singleAssoc.txt.gz", old_sa)
w("in/fxv4-test.pooledTests.txt.gz", ["Region\tRows\tRefused\tRule", "R1\t3\t0\tbonferroni", "R3\t1\t0\tcauchy"])
for tag, (res, sa) in plants.items():
    os.makedirs(f"in/{tag}", exist_ok=True); w(f"in/{tag}/fxv5-test.txt.gz", res); w(f"in/{tag}/fxv5-test.singleAssoc.txt.gz", sa)
    gzip.open(f"in/{tag}/fxv5-test.pooledTests.txt.gz", "wt").write("# placeholder\n")
w("in/ew/fxv5-test.exactByWeight.txt.gz", ["Kind\tRegion\tGroup\tmax_MAF\tWeight\tMAC\tn_eff\tn_material", "burden\tR2\tgrp\t0.001\tW1\t7\t1.0000012\t1"])
PH = "Model\tKind\tRegion\tGroup\tmax_MAF\tWeight\tMarkerID\tn_dropped\tdropped_mass\tguard"
w("in/pn/fxv5-test.spaPinned.txt.gz",     [PH, "1\tburden\tR2\tgrp\t0.001\tW1\t\t41\t4.1e-11\tapplied", "1\tburden\tR3\tgrp\t0.001\tW1\t\t41\t4.1e-11\tapplied"])
w("in/pnref/fxv5-test.spaPinned.txt.gz",  [PH, "1\tburden\tR2\tgrp\t0.001\tW1\t\t0\t1.1e-06\trefused"])
w("in/pnreal/fxv5-test.spaPinned.txt.gz", [PH, "1\tsingle_in_group\tR1\t\t\t\tchr1:300:A:G\t22\t2.2e-11\tapplied"])
SG = "Model\tKind\tRegion\tGroup\tmax_MAF\tWeight\tgate\tz\tp_spa\tp_exact\treported"
w("in/sg/fxv5-test.stretchGate.txt.gz",    [SG, "1\tburden\tR2\tgrp\t0.001\tW1\t1.61\t3.8\t0.3\t0.31\texact", "1\tburden\tR3\tgrp\t0.001\tW1\t0.84\t3.6\t0.5\t0.51\texact"])
w("in/sgspa/fxv5-test.stretchGate.txt.gz", [SG, "1\tburden\tR2\tgrp\t0.001\tW1\t0.12\t2.3\t0.3\tNA\tsaddlepoint"])
SG9 = "Model\tKind\tRegion\tGroup\tmax_MAF\tWeight\tgate\tz\tp_spa"
w("in/sghdr/fxv5-test.stretchGate.txt.gz", [SG9, "1\tburden\tR1\tgrp\t0.001\tW3\t0.15\t2.4\t0.4",
                                            "1\tburden\tR2\tgrp\t0.001\tW1\t1.61\t3.8\t0.3\t0.31\texact"])
w("in/pnur/fxv5-test.spaPinned.txt.gz",   [PH, "1\tultra_rare\tR6\tgrp\t0.001\tAM_cbrt0.2_unw\tR6:grp:0.001:AM_cbrt0.2_unw\t26\t2.6e-11\tapplied"])
# a LONG list (20,000 rows of regions that are not in the table, ~1 MB): a `| head -1` on it dies of SIGPIPE under pipefail
w("in/allok/fxv5-test.exactByWeight.txt.gz", ["Kind\tRegion\tGroup\tmax_MAF\tWeight\tMAC\tn_eff\tn_material", "burden\tR2\tgrp\t0.001\tW1\t7\t1.0000012\t1"]
  + [f"burden\tZ{i:05d}\tgrp\t0.001\tAM_sbeta20_unw\t30\t1.0000001\t1" for i in range(20000)])
w("in/nosingle/fxv5-test.exactByWeight.txt.gz", ["Kind\tRegion\tGroup\tmax_MAF\tWeight\tMAC\tn_eff\tn_material", "burden\tR2\tgrp\t0.001\tW1\t7\t1.0000012\t1"])
for tag in ("allok", "allbad", "nosingle"):   # a complete rerun discloses the same pooled rows as the old run
    w(f"in/{tag}/fxv5-test.pooledTests.txt.gz", ["Region\tRows\tRefused\tRule", "R1\t3\t0\tbonferroni", "R3\t1\t0\tcauchy"])
w("in/regions.txt", ["R1", "R2", "R4", "R5", "R6", "R7"], gz=False)
w("in/failures.tsv", ["file\tMarkerID\tRegion\tCHR\tp.value\tp.value.NA\tIs.SPA",
                      "fxv4-test\tR1:grp:0.001:W1\tR1\t1\t1e-100\t1e-100\tfalse",
                      "fxv4-test\tchr1:100:A:G\t-\t1\t1e-90\t1e-90\tfalse"], gz=False)
# stage-checker fixtures
os.makedirs("stage/null", exist_ok=True); os.makedirs("stage/pos", exist_ok=True); os.makedirs("stage/chr", exist_ok=True)
hdr = "Kind\tRegion\tGroup\tmax_MAF\tWeight\tMarkerID\tMAC\tp_normal_not_reported\tp_exact_reported\tNote"
w("stage/null/a.spaFallbacks.txt.gz", [hdr])
w("stage/pos/b.spaFallbacks.txt.gz", [hdr, "region\tG1\tpLoF\t0.001\tunw_unw\tNA\t9\t1e-120\t0.0021\tspa: root(q) not found [root1=inf/1000 root2=-3/5 qsaddle=n/a bound=-1]",
                                          "region\tG3\tpLoF\t0.001\tunw_unw\tNA\t44\t1e-200\tNA\tspa: saddle(q) invalid [qsaddle=INVALID]; exact test not run or failed: above the ceiling"])
w("stage/chr/c.spaFallbacks.txt.gz", [hdr]); w("stage/chr/c.chrStatus.txt.gz", ["requested\t1 2 3", "failed\t5", "refused_by_gate\tnone", "skipped_no_region\tnone"])
w("stage/null/a.chrStatus.txt.gz", ["requested\t1 2 3", "failed\tnone", "refused_by_gate\tnone", "skipped_no_region\tnone"])
print("fixtures written")
PY
runall(){ local tag="$1" newtag="$2"; rm -rf "work/${tag}" "out/${tag}"; mkdir -p "work/${tag}" "out/${tag}"
  (cd "work/${tag}" && OLD_PREFIX="${T}/in/fxv4-test" NEW_PREFIX="${T}/in/${newtag}/fxv5-test" REGIONS=ALL FAILURES=NONE OUT_DIR="${T}/out/${tag}" bash "${HERE}/aou_spa_splice.sh" > "${T}/${tag}.log" 2>&1); echo $?; }
expectall(){ local tag="$1" newtag="$2" want="$3" pat="$4"; local got; got=$(runall "${tag}" "${newtag}")
  if [[ ${got} != "${want}" ]] || ! grep -q -- "${pat}" "${tag}.log"; then echo "FAIL ${tag}: exit ${got} (want ${want}), pattern '${pat}' $(grep -q -- "${pat}" "${tag}.log" && echo found || echo MISSING)"; fail=1; else echo "ok   ${tag}: exit ${got}, '${pat}'"; fi; }
run(){ local tag="$1"; rm -rf "work/${tag}" "out/${tag}"; mkdir -p "work/${tag}" "out/${tag}"
  (cd "work/${tag}" && OLD_PREFIX="${T}/in/fxv4-test" NEW_PREFIX="${T}/in/${tag}/fxv5-test" REGIONS="${T}/in/regions.txt" FAILURES="${T}/in/failures.tsv" OUT_DIR="${T}/out/${tag}" bash "${HERE}/aou_spa_splice.sh" > "${T}/${tag}.log" 2>&1); echo $?; }
fail=0
expect(){ local tag="$1" want="$2" pat="$3"; local got; got=$(run "${tag}")
  if [[ ${got} != "${want}" ]] || ! grep -q -- "${pat}" "${tag}.log"; then echo "FAIL ${tag}: exit ${got} (want ${want}), pattern '${pat}' $(grep -q -- "${pat}" "${tag}.log" && echo found || echo MISSING)"; fail=1; else echo "ok   ${tag}: exit ${got}, '${pat}'"; fi; }
expect ok   0 "unobservable burden cells (Number_rare > 0) accepted at region level: 2 (p rose in 2, fell in 0)"
expect ew   0 "cells routed to the exact test by weight (LEDGER #172): 1"
grep -q "ROWS .pooledTests.txt.gz: R1 old 1 new 0 (a region with nothing to disclose" ok.log && echo "ok   row counts per region: the ROWS line of the ok case" || { echo "FAIL the ROWS line is missing from ok.log"; fail=1; }
grep -c "^ROWS " allok.log >/dev/null 2>&1 || true
grep -q "pseudo-marker rows explained by their own exact-test eligibility (LEDGER #166, count <= 4): 4" ok.log || { echo "FAIL ok: pseudo-marker count"; fail=1; }
grep -q "real-marker rows: 4 rerun marker(s), 4 old row(s) replaced by 4 new, 2 old cop(y/ies) kept" <(tr -s ' ' < ok.log) || { echo "FAIL ok: real-marker accounting"; fail=1; }
[[ $(gzip -dc out/ok/fxv4-test.singleAssoc.txt.gz | grep -c "chr1:500") -eq 3 ]] || { echo "FAIL ok: shared marker copies"; fail=1; }
expect real 1 "ERROR .singleAssoc.txt.gz"
expect out  1 "new rows outside the region list"
expect obs  1 "ERROR .txt.gz"
expect dup  1 "ERROR .singleAssoc.txt.gz"
expect ur   0 "pseudo-marker rows accepted on STATISTIC IDENTITY ALONE"
expect twin 0 "pseudo-marker rows accepted on STATISTIC IDENTITY ALONE"
grep -q "columns that differ old -> new" ur.log && { echo "FAIL ur: a row accepted on statistic identity must not be reported as unexplained"; fail=1; } || echo "ok   ur: accepted without an unexplained-change diagnostic"
[[ -s out/ur/fxv4-test.spliceStatId.txt ]] && echo "ok   the statistic-identity list travels with the spliced file" || { echo "FAIL ur: .spliceStatId.txt missing"; fail=1; }
grep -q "rows accepted on statistic identity alone (CONTROL 1c, the weaker tier): 1" out/ur/fxv4-test.splice.txt || { echo "FAIL ur: manifest count"; fail=1; }
expect sidvar  1 "ERROR .singleAssoc.txt.gz"
expect sidbeta 1 "ERROR .singleAssoc.txt.gz"
expect sidac   1 "ERROR .singleAssoc.txt.gz"
grep -q "columns that differ old -> new: var" sidvar.log || { echo "FAIL sidvar: the refusal must name the differing column"; fail=1; }
expect pn     0 "cells in which pinned samples were left out of the saddlepoint (LEDGER #174): 1"
grep -q ".spaPinned.txt.gz: 1 guarded test(s) in the spliced regions, carried over" pn.log && echo "ok   the pinned sidecar is carried over restricted to the spliced regions (R3 is not one)" || { echo "FAIL pn: sidecar carry-over"; fail=1; }
expect pnref  1 "ERROR .txt.gz"
expect sg     0 "cells reported by the exact convolution at the stretch gate (LEDGER #176): 1"
grep -q ".stretchGate.txt.gz: 1 gated cell(s) in the spliced regions, carried over" sg.log && echo "ok   the stretch-gate sidecar is carried over restricted to the spliced regions (R3 is not one)" || { echo "FAIL sg: sidecar not carried"; fail=1; }
expect sgspa  1 "ERROR .txt.gz"
expect sghdr  0 "cells reported by the exact convolution at the stretch gate (LEDGER #176): 1"
expect pnreal 0 "real markers in which pinned samples were left out of the saddlepoint (LEDGER #174): 1"
expect pnur   0 "cells in which pinned samples were left out of the saddlepoint (LEDGER #174): 1"
expectall allok  allok  0 "cells routed to the exact test by weight (LEDGER #172): 1"
grep -q ".exactByWeight.txt.gz: 1 routed cell(s) in the spliced regions, carried over" allok.log && [[ $(gzip -dc out/allok/fxv4-test.exactByWeight.txt.gz | wc -l) -eq 2 ]] && echo "ok   a 20,000-row routed list is carried over restricted to the table's regions" || { echo "FAIL the routed list carry-over"; fail=1; }
expectall allbad allbad 1 "ERROR .txt.gz"
expectall nosingle nosingle 0 "NOTE .singleAssoc.txt.gz: the rerun had SINGLE_IN_GROUP=FALSE (placeholder); not spliced, not carried over"
[[ ! -e out/nosingle/fxv4-test.singleAssoc.txt.gz ]] && echo "ok   nosingle: no single-variant file written for the splice" || { echo "FAIL nosingle: a single-variant file was written from a placeholder"; fail=1; }
# stage checker
bash "${HERE}/aou_spa_stagecheck.sh" stage/null > sc_null.log 2>&1 && grep -q "STAGE OK" sc_null.log && echo "ok   stagecheck null" || { echo "FAIL stagecheck null"; fail=1; }
if bash "${HERE}/aou_spa_stagecheck.sh" stage/pos > sc_pos.log 2>&1; then echo "FAIL stagecheck pos: exit 0"; fail=1; else grep -q "NA 1" sc_pos.log && echo "ok   stagecheck NA flagged" || { echo "FAIL stagecheck pos"; fail=1; }; fi
if bash "${HERE}/aou_spa_stagecheck.sh" stage/chr > sc_chr.log 2>&1; then echo "FAIL stagecheck chr: exit 0"; fail=1; else grep -q "CHROMOSOME FAILED" sc_chr.log && echo "ok   stagecheck lost chromosome flagged" || { echo "FAIL stagecheck chr"; fail=1; }; fi
if (( fail )); then echo "SELFTEST FAILED"; exit 1; fi
echo "selftest OK: 11 positive controls (two via the exactByWeight sidecar, three via the spaPinned sidecar, two via the stretchGate sidecar (one under a header that lacks the column), two on statistic identity, two complete reruns -- one with the single-variant file switched off), 10 planted refusals (one a guard that refused, three moving a statistic the exact test never touches), 3 stage-checker cases"
