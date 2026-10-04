#!/bin/bash
# Splice a filtered FlexRV rerun (LEDGER #164) over the original delivered files.
# Runs INSIDE a dsub job (gcr.io/google.com/cloudsdktool/google-cloud-cli:slim);
# streams with gcloud storage cat. Env:
#   OLD_PREFIX  gs://.../<pre>          (v4: .txt.gz, .singleAssoc.txt.gz, .pooledTests.txt.gz)
#   NEW_PREFIX  gs://.../<pre>          (v5 rerun of the listed regions only)
#   REGIONS     gs://.../rerun/<pre>.regions.txt, or ALL for a COMPLETE rerun: every
#               region of the new results file (LEDGER #172's unfiltered rerun;
#               derived inside the job, so the notebook reads no results file)
#   FAILURES    gs://.../failures.tsv   (the scan's; rows of this file's failures are
#               used), or NONE when OLD_PREFIX already carries the earlier fixes
#   OUT_DIR     local dir dsub delocalizes (--output-recursive)
set -euo pipefail
: "${OLD_PREFIX:?}"; : "${NEW_PREFIX:?}"; : "${REGIONS:?}"; : "${FAILURES:?}"; : "${OUT_DIR:?}"
mkdir -p "${OUT_DIR}"; pre=$(basename "${OLD_PREFIX}")
cat_any(){  # gs:// or a local path, gzipped or not
  if [[ $1 == gs://* ]]; then src=(gcloud storage cat "$1"); else src=(cat "$1"); fi
  if [[ $1 == *.gz ]]; then "${src[@]}" | gzip -dc; else "${src[@]}"; fi; }
exists(){ if [[ $1 == gs://* ]]; then gcloud storage ls "$1" >/dev/null 2>&1; else [[ -s $1 ]]; fi; }
if [[ ${REGIONS} == "ALL" ]]; then
  cat_any "${NEW_PREFIX}.txt.gz" | awk -F'\t' 'NR>1 {print $1}' | sed '/^$/d' | sort -u > regions.txt
  echo "REGIONS=ALL: every region of the new results file"
else
  cat_any "${REGIONS}" | sed '/^$/d' | sort -u > regions.txt
fi
if [[ ${FAILURES} == "NONE" ]]; then
  printf 'file\tMarkerID\tRegion\tCHR\tp.value\tp.value.NA\tIs.SPA\n' > failures.none.tsv
  FAILURES="${PWD}/failures.none.tsv"
  echo "FAILURES=NONE: no recorded failure is expected to explain a change"
fi
cat_any "${FAILURES}" | awk -F'\t' -v f="${pre}" 'NR>1 && $1==f && $3!="-" {print $3}' | sort -u > failregions.txt
echo "regions to splice: $(wc -l < regions.txt); with failed rows: $(wc -l < failregions.txt)"

# region of a row: results/pooled = column Region; singleAssoc = MarkerID prefix (UR rows)
splice_one(){  # suffix  region-column-rule (awk expr on $0 giving the region)
  local sfx="$1" rule="$2"
  local old="${OLD_PREFIX}${sfx}" new="${NEW_PREFIX}${sfx}"   # a separate `local`: bash expands a
  local out="${OUT_DIR}/${pre}${sfx%.gz}"                     # local's words before it assigns them
  cat_any "${old}" > old.tsv; cat_any "${new}" > new.tsv
  # A sidecar the run had nothing to write is delivered as a one-line '#'
  # placeholder (the wrapper's doing: dsub fails delocalization on a declared
  # output that does not exist). Treat a placeholder as an EMPTY table with the
  # other side's header; both placeholders means nothing to splice.
  local oph=0 nph=0
  [[ $(head -c 1 old.tsv) == "#" ]] && oph=1
  [[ $(head -c 1 new.tsv) == "#" ]] && nph=1
  if (( oph && nph )); then
    cp old.tsv "${out}"; gzip -f "${out}"; echo "${sfx}: placeholder on both sides, carried over"; return 0
  elif (( oph )); then head -1 new.tsv > old.tsv; echo "${sfx}: old is a placeholder (no region had anything to disclose); treating as empty"
  elif (( nph )); then head -1 old.tsv > new.tsv; echo "${sfx}: new is a placeholder; treating as empty"
  fi
  head -1 old.tsv > hdr.old; head -1 new.tsv > hdr.new
  cmp -s hdr.old hdr.new || { echo "ERROR ${sfx}: header differs between old and new" >&2; exit 1; }
  # The singleAssoc holds TWO kinds of row (found on the eas/afr stage,
  # 2026-09-27). A pseudo-marker row's MarkerID is Region:Group:maxMAF:Weight
  # and names its region. A REAL marker -- rare but not ultra-rare, MAC above
  # the collapse cutoff, which exists whenever 2*N*maxMAF > 10, i.e. N above
  # ~5,000 at max MAF 0.001: eas and afr, never sas -- is tested singly under
  # its own chr:pos:ref:alt ID, which names no region. Real-marker rows are
  # keyed by MarkerID: the rerun's rows replace the old rows with the same ID
  # (the rerun processed only the listed regions, so those ARE the rerun
  # regions' markers), and a changed one must be a recorded real-marker
  # failure or a marker at MAC <= ER_MAC. A row is a real-marker row when the
  # MarkerID's second ':'-field is a number (a position).
  local single=0; [[ ${sfx} == ".singleAssoc.txt.gz" ]] && single=1
  local isreal='(split($3, a_, ":") >= 2 && a_[2] ~ /^[0-9]+$/)'
  # newmarkers.txt: MarkerID<TAB>copies in new. A real marker sits in every
  # region that spans it (the protocadherin clusters put one variant in a
  # dozen), and its row is written once per region, so the old file can hold
  # k identical copies where the rerun of one region writes one: MULTISET
  # semantics -- replace as many old copies as new has, keep the rest for the
  # regions not rerun (afr/VTE chr5:140.98 Mb, 56 surplus copies refused
  # before this).
  : > newmarkers.txt
  if (( single )); then awk -F'\t' 'FNR>1 && '"${isreal}"' {c[$3]++} END{for(m in c) print m"\t"c[m]}' new.tsv | sort > newmarkers.txt; fi
  # new pseudo-marker/region rows, keyed by region (only the listed regions are expected in new)
  awk -F'\t' -v rule="${rule}" -v single="${single}" 'NR==FNR{keep[$1]=1; next} FNR==1{next} single && '"${isreal}"' {next} {reg=('"${rule}"'); if(!(reg in keep)){bad++} } END{if(bad){print "ERROR: " bad " new rows outside the region list" > "/dev/stderr"; exit 1}}' regions.txt new.tsv
  # CONTROL 1 (row-level, since LEDGER #166): a rerun changes exactly two
  # kinds of row -- a cell whose saddlepoint failed under v4 (in failures.tsv)
  # and a cell at MAC <= --max_MAC_for_ER (ER_MAC, default 4), which now takes
  # the exact test whether or not it failed. Any OTHER changed row is a
  # refusal: it would mean the rerun differs from v4 for a reason we have not
  # named. Cells are keyed Region|Weight (results/pooled: the row's own
  # columns; singleAssoc pseudo-marker rows: the MarkerID's first and last
  # fields); a real-marker row is keyed by its MarkerID. MAC is read from the
  # OLD file: per set from the results, per marker from AC_Allele2 folded by
  # AF_Allele2 in the singleAssoc.
  ERMAC="${ER_MAC:-4}"
  local nri=0
  if [[ ${sfx} == ".txt.gz" ]]; then
    awk -F'\t' 'FNR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["Group"]!="Cauchy"{print $h["Region"]"|"$h["Weight"]"\t"$h["MAC"]}' old.tsv > cellmac.tsv
    nri=$(head -1 old.tsv | tr '\t' '\n' | awk '$0=="Number_rare"{print NR}'); nri=${nri:-0}
    pci=$(head -1 old.tsv | tr '\t' '\n' | awk '$0=="Pvalue_Burden"{print NR}'); pci=${pci:-5}
    : > okregions.txt
  fi
  : > markermac.tsv; : > urown.tsv
  if (( single )); then
    awk -F'\t' 'FNR==1{for(i=1;i<=NF;i++)h[$i]=i; next} '"${isreal}"' {ac=$h["AC_Allele2"]+0; af=$h["AF_Allele2"]+0; mac=ac; if(af>0.5 && af<1) mac=ac*(1-af)/af; print $3"\t"mac}' old.tsv > markermac.tsv
    # A pseudo-marker's single test takes the exact test when ITS OWN unweighted
    # allele count is <= ER_MAC (LEDGER #166 gave that to every set; v4 gave it
    # to unw_unw alone), and that count is not the set's MAC in the results
    # (which adds the rare, non-collapsed markers) nor the row's AC_Allele2 for
    # a weighted set (a weighted sum). The exact handle: a set with 0/1 weights
    # (unw|lof score x unw|MACk MAF) reports its own count (LEDGER #133), and
    # every set collapses the same markers as -- or a subset of -- the 0/1 set
    # with the same zero pattern (lof x mafbeta keeps the LoF markers, as
    # lof x unw does; cbrt x MAC25 keeps the MAC <= 25 markers, as unw x MAC25
    # does), so that twin's count bounds it from above. Found on eas/Asth and afr/BenIntNeo (a lof_MAC100
    # pseudo-marker whose own count is at or below the cutoff, inside a set
    # whose MAC is above it).
    awk -F'\t' 'FNR==1{for(i=1;i<=NF;i++)h[$i]=i; next} !('"${isreal}"') {n=split($3,a,":"); if(n<2) next; ac=$h["AC_Allele2"]+0; af=$h["AF_Allele2"]+0; mac=ac; if(af>0.5 && af<1) mac=ac*(1-af)/af
        if(a[n] ~ /_(unw|lof)_(unw|MAC[0-9]+)$/) print a[1]"|"a[n]"\t"mac > "urown.tsv" }' old.tsv
  fi
  cat_any "${FAILURES}" | awk -F'\t' -v f="${pre}" 'NR>1 && $1==f {n=split($2,a,":"); if(!(n>=2 && a[2] ~ /^[0-9]+$/)) print a[1]"|"a[n]}' | sort -u > failcells.txt
  cat_any "${FAILURES}" | awk -F'\t' -v f="${pre}" 'NR>1 && $1==f {n=split($2,a,":"); if(n>=2 && a[2] ~ /^[0-9]+$/) print $2}' | sort -u > failmarkers.txt
  # LEDGER #172 (2026-09-29): a cell the rerun routed to the exact test because
  # one atom carries its weight is named in the rerun's own
  # `.exactByWeight.txt.gz` (Kind Region Group max_MAF Weight MAC n_eff
  # n_material); its change is explained by that listing, Region|Weight.
  : > ewcells.txt
  if exists "${NEW_PREFIX}.exactByWeight.txt.gz"; then
    cat_any "${NEW_PREFIX}.exactByWeight.txt.gz" | awk -F'\t' 'NR>1 {print $2"|"$5}' | sort -u > ewcells.txt
  fi
  # LEDGER #176 (FIXED 2026-10-03): a cell the stretch gate sent to the exact
  # convolution is named in the rerun's own `.stretchGate.txt.gz` (Model Kind
  # Region Group max_MAF Weight gate z p_spa p_exact reported ...). Only a row
  # whose `reported` reads `exact` explains a change: `saddlepoint` means the
  # cell was evaluated and listed (at or above the floor) but the p is the one
  # it always was. Region|Weight, as the routed cells.
  #
  # The sidecar the wrapper delivers is the per-chromosome files concatenated
  # under the FIRST file's header, and the image that ran the v8 pass
  # (1.5.2-dev-152ffd8c) writes the `p_exact` and `reported` columns only when
  # a chromosome had a gated cell. So the header can lack `reported` while
  # later rows carry it in field 11 -- and an awk lookup of a missing header
  # name is `$0`, which matches nothing, which refused 25 of afr's 124 checks
  # on 2026-10-04. Take the column by name when the header has it and by its
  # position when it does not; a 9-field row has no field 11 and is never
  # accepted, which is right (it was not gated).
  : > sgcells.txt
  if exists "${NEW_PREFIX}.stretchGate.txt.gz"; then
    cat_any "${NEW_PREFIX}.stretchGate.txt.gz" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; c=("reported" in h)?h["reported"]:11; next} (NF>=c && $c=="exact") {print $3"|"$6}' | sort -u > sgcells.txt
  fi
  # LEDGER #174 (2026-09-29): a test in which pinned samples were left out of
  # the saddlepoint's vector is named in the rerun's own `.spaPinned.txt.gz`
  # (Model Kind Region Group max_MAF Weight MarkerID n_dropped dropped_mass
  # guard). Only a row whose guard reads `applied` explains a change: `refused`
  # means the mass was over the budget, nothing was left out and the test ran
  # as before. A cell is keyed Region|Weight (kinds burden, ultra_rare and
  # region_rescale); a real marker by its MarkerID (single_in_group, marker).
  : > pncells.txt; : > pnmarkers.txt
  if exists "${NEW_PREFIX}.spaPinned.txt.gz"; then
    cat_any "${NEW_PREFIX}.spaPinned.txt.gz" | awk -F'\t' 'NR>1 && $10=="applied" {
        n=split($7, a, ":")
        if(($2=="single_in_group" || $2=="marker") && n>=2 && a[2] ~ /^[0-9]+$/) print $7 > "pnmarkers.raw"
        else if($3!="") print $3"|"$6 > "pncells.raw" }'
    [[ -f pncells.raw ]] && { sort -u pncells.raw > pncells.txt; rm -f pncells.raw; }
    [[ -f pnmarkers.raw ]] && { sort -u pnmarkers.raw > pnmarkers.txt; rm -f pnmarkers.raw; }
  fi
  # changed rows = old rows the rerun replaces (listed regions' rows; for the
  # singleAssoc also the rerun markers' rows) that do not appear verbatim in new, and vice versa
  awk -F'\t' -v rule="${rule}" -v single="${single}" 'FILENAME=="regions.txt"{keep[$1]=1; next} FILENAME=="newmarkers.txt"{nm[$1]=$2; next} FNR==1{next} { if(single && '"${isreal}"'){ if($3 in nm){ seen[$3]++; if(seen[$3] <= nm[$3]) print } next } reg=('"${rule}"'); if(reg in keep) print }' regions.txt newmarkers.txt old.tsv | sort > o.all
  awk -F'\t' 'FNR>1' new.tsv | sort > n.all
  comm -3 o.all n.all | sed 's/^\t//' > changed.rows
  # CONTROL 1c (2026-09-30): STATISTIC IDENTITY -- the LAST tier, a deliberately
  # WEAKER claim than every rule above, and singleAssoc PSEUDO-MARKER rows only.
  #
  # Why it exists. #166 eligibility is `MAC_nw <= ER_MAC`, the set's UNWEIGHTED
  # collapse count (src/Main.cpp:4132-4134). That number is in NO delivered
  # file when the cell has non-collapsed rare markers: the results `MAC` column
  # is MAC_nw PLUS the kept rare markers' MACs, and `.markerList.txt` is per
  # CELL, not per SET. The twin rule above bounds it from above and the bound is
  # loose exactly when the twin keeps markers the set drops -- a gene whose rare
  # variants are largely unscored by AlphaMissense, where `unw` keeps every
  # allele and a nonlinear transform keeps only the scored ones. That refused
  # four prefixes of the 2026-09-30 rollout (amr Pancreat/Stroke, afr Asth/VTE).
  #
  # DO NOT "FIX" THAT BY READING THE ROW'S OWN AC_Allele2. A pseudo-marker's
  # AC_Allele2 is the WEIGHTED collapse (`paritySumDense(genoURVec)`,
  # src/Main.cpp:3888), which #118 deliberately left weighted -- that is why
  # every `.singleAssoc.txt` was byte-identical through that fix. MEASURED on
  # the committed golden `plink_group_flexrv_skato.txt`: GENE1's four sets have
  # identical per-set unweighted composition (MAC 183, Number_rare 6,
  # Number_ultra_rare 17 -- all per-set since #118) and report AC_Allele2 54,
  # 50.47, 47.85 and 44.4398. Comparing that against an allele-count cutoff is
  # LEDGER #134's defect (the SCALE of the weights choosing the statistical
  # treatment) moved into the control, and it is ANTICONSERVATIVE: every
  # down-scaling transform reads low, so it would accept rows that are not
  # #166-eligible. The 0/1 twin reads an integer, which is why the twin works.
  #
  # What this tier asserts instead, MEASURED (2026-09-30) on
  # input/group_new_snpid_withWeights_atom.txt:
  #   POSITIVE, #166 arm  (--max_MAC_for_ER 0 -> 50):     8 rows move, columns {SE, p.value}
  #   POSITIVE, #172 arm  (SAIGE_EXACT_BY_WEIGHT 0 -> 1): 10 rows move, columns {SE, p.value}
  #   NEGATIVE, variance ratio perturbed by 0.1%:         14 rows move, columns {BETA, SE, p.value, p.value.NA, var}
  # So a routing change moves SE and p.value and NOTHING else: BETA and Tstat
  # are the score statistic, `var` its variance, AC_Allele2/AF_Allele2 the
  # genotype summary, and the exact test re-derives SE from the new p (the
  # eas/T2Diab reading, "same BETA, SE re-derived from the new p"). A change in
  # the DATA or the MODEL moves `var` -- note Tstat alone does NOT discriminate,
  # since the score does not depend on the variance ratio.
  #
  # THE LIMIT, stated so it is not mistaken for the rules above: this is a
  # FINGERPRINT, not a REASON. It says the statistic is untouched and only the
  # p-value method moved; it does NOT say the method chosen was the right one,
  # so a cell that switched methods WRONGLY (LEDGER #176's shape) would pass.
  # Every row it accepts is counted, named on stderr and written to
  # `<pre>.spliceStatId.txt` beside the spliced file. The exact fix is a sidecar
  # from the RUN naming every #166-routed row with its MAC_nw, as
  # `.exactByWeight.txt` and `.spaPinned.txt` already do for #172 and #174.
  : > rowdiff.tsv; : > statid.txt
  if (( single )); then
    awk -F'\t' '
      FILENAME=="hdr.old"{ for(i=1;i<=NF;i++) nm[i]=$i; nf=NF; next }
      FILENAME=="o.all"{ if(!('"${isreal}"')) o[$3]=$0; next }
      { if('"${isreal}"') next; if(!($3 in o)) next
        split(o[$3], a, "\t"); d=""
        for(i=1;i<=nf;i++) if(a[i] != $i) d = (d=="" ? nm[i] : d","nm[i])
        if(d!="") print $3"\t"d }' hdr.old o.all n.all > rowdiff.tsv
    awk -F'\t' '{ n=split($2,c,","); keep=1; for(i=1;i<=n;i++) if(c[i]!="SE" && c[i]!="p.value") keep=0
                  if(keep) print $1 }' rowdiff.tsv | sort -u > statid.txt
  fi
  # A burden cell with rare, NON-collapsed markers (Number_rare > 0) is a test
  # the scan could not observe: its vector is w*(real markers) + pseudo-marker
  # and its saddlepoint can fail while every component's own single test
  # converged, and v4 recorded nothing (pre-#164 the normal p was reported
  # silently). Found on eas/T2Diab: a MAC-47 sbeta20 cell at 1.6e-7 in v4 and
  # 0.068 on the fixed image, same BETA, SE re-derived from the new p. Such a
  # cell's change is accepted at REGION level (the region is in the rerun list
  # for a recorded failure or a tail p, which is what the tail arm exists
  # for), counted, and reported with its direction; the exact cell-level rule
  # stays for every observable cell (Number_rare == 0, real markers, pseudo-
  # marker rows). Pass 1 judges the cells and records the regions with an
  # explained change; pass 2 judges the Cauchy/pooled rows against them.
  awk -F'\t' -v sfx="${sfx}" -v single="${single}" -v nri="${nri}" 'FILENAME=="failcells.txt"{fc[$1]=1; n=split($1,b,"|"); fr[b[1]]=1; next}
        FILENAME=="failmarkers.txt"{fm[$1]=1; next}
        FILENAME=="cellmac.tsv"{mac[$1]=$2; if($2+0 <= '"${ERMAC}"'){ split($1,b,"|"); mr[b[1]]=1 } next}
        FILENAME=="markermac.tsv"{mmac[$1]=$2; next}
        FILENAME=="urown.tsv"{uown[$1]=$2; next}
        FILENAME=="okregions.txt"{okr[$1]=1; next}
        FILENAME=="ewcells.txt"{ew[$1]=1; next}
        FILENAME=="sgcells.txt"{sg[$1]=1; next}
        FILENAME=="pncells.txt"{pn[$1]=1; next}
        FILENAME=="pnmarkers.txt"{pnm[$1]=1; next}
        FILENAME=="statid.txt"{sid[$1]=1; next}
        FILENAME=="rowdiff.tsv"{rd[$1]=$2; next}
        FNR==1 && FILENAME=="changed.rows"{pass++}
        pass==1 {
          if(sfx==".pooledTests.txt.gz" || $2=="Cauchy") next
          ok=0; unobs=0
          if(single && '"${isreal}"'){ key=$3; reg=""; if(key in fm) ok=1; if((key in mmac) && mmac[key]+0 <= '"${ERMAC}"') ok=1
                                       if(!ok && (key in pnm)){ ok=1; pnmark[key]=1 } }
          else if(single){ id=$3; n=split(id,a,":"); key=a[1]"|"a[n]; reg=a[1]; if(key in fc) ok=1; if((key in mac) && mac[key]+0 <= '"${ERMAC}"') ok=1
                           if(!ok && (key in ew)){ ok=1; ewcell[key]=1 }
                           if(!ok && (key in sg)){ ok=1; sgcell[key]=1 }
                           if(!ok && (key in pn)){ ok=1; pncell[key]=1 }
                           if(!ok){ # the 0/1-weight twin with the same zero pattern collapses the same markers (or a superset): lof stays lof, any other score -> unw; MACk stays, any other MAF -> unw
                             w=a[n]; m=split(w,f,"_"); if(m>=3){ sc=f[m-1]; mf=f[m]; pre=f[1]; for(i=2;i<=m-2;i++) pre=pre"_"f[i]
                               twin=a[1]"|"pre"_"((sc=="lof")?"lof":"unw")"_"((mf ~ /^MAC[0-9]+$/)?mf:"unw")
                               if((twin in uown) && uown[twin]+0 <= '"${ERMAC}"'){ ok=1; urer[key]=1 } } }
                           if(!ok && ($3 in sid)){ ok=1; sidrow[$3]=1 } }
          else { key=$1"|"$4; reg=$1; if(key in fc) ok=1; if((key in mac) && mac[key]+0 <= '"${ERMAC}"') ok=1
                 if(!ok && (key in ew)){ ok=1; ewcell[key]=1 }
                 if(!ok && (key in sg)){ ok=1; sgcell[key]=1 }
                 if(!ok && (key in pn)){ ok=1; pncell[key]=1 }
                 if(!ok && nri>0 && $nri+0 > 0){ ok=1; unobs=1; ucell[key]=1 } }
          if(ok){ if(reg!="") okr[reg]=1 } else { bad++
            if(bad<=5){ print "  unexplained change: " substr($0,1,120) > "/dev/stderr"
                        if(single && ($3 in rd)) print "    columns that differ old -> new: " rd[$3] > "/dev/stderr" } }
          next }
        pass==2 {
          if(!(sfx==".pooledTests.txt.gz" || $2=="Cauchy")) next
          reg=$1; ok=0; if(reg in okr) ok=1; if(reg in fr) ok=1; if(reg in mr) ok=1
          if(!ok){ bad++; if(bad<=5) print "  unexplained change: " substr($0,1,120) > "/dev/stderr" } }
        END{ for(k in ucell) print k > "ukeys.tmp"
             nu=0; for(k in urer) nu++; if(nu>0) print "  pseudo-marker rows explained by their own exact-test eligibility (LEDGER #166, count <= '"${ERMAC}"'): " nu > "/dev/stderr"
             ewn=0; for(k in ewcell) ewn++; if(ewn>0) print "  cells routed to the exact test by weight (LEDGER #172): " ewn > "/dev/stderr"
             sgn=0; for(k in sgcell) sgn++; if(sgn>0) print "  cells reported by the exact convolution at the stretch gate (LEDGER #176): " sgn > "/dev/stderr"
             pnn=0; for(k in pncell) pnn++; if(pnn>0) print "  cells in which pinned samples were left out of the saddlepoint (LEDGER #174): " pnn > "/dev/stderr"
             pmn=0; for(k in pnmark) pmn++; if(pmn>0) print "  real markers in which pinned samples were left out of the saddlepoint (LEDGER #174): " pmn > "/dev/stderr"
             sdn=0; for(k in sidrow){ sdn++; print k > "statid.accepted" }
             if(sdn>0) print "  pseudo-marker rows accepted on STATISTIC IDENTITY ALONE (only SE and p.value moved; BETA, Tstat, var and AC_Allele2 byte-identical -- the cell was NOT named, see CONTROL 1c): " sdn > "/dev/stderr"
             for(r in okr) print r > "okregions.tmp"
             print bad+0 }' failcells.txt failmarkers.txt cellmac.tsv markermac.tsv urown.tsv okregions.txt ewcells.txt sgcells.txt pncells.txt pnmarkers.txt statid.txt rowdiff.tsv changed.rows changed.rows > nbad.txt
  nbad=$(cat nbad.txt); [[ -f okregions.tmp ]] && { cat okregions.tmp >> okregions.txt; rm -f okregions.tmp; sort -u -o okregions.txt okregions.txt; }
  if [[ -f statid.accepted ]]; then   # provenance travels with the file, not just the log
    sort -u statid.accepted >> "${OUT_DIR}/${pre}.spliceStatId.txt"; rm -f statid.accepted
  fi
  if [[ -s ukeys.tmp ]]; then   # direction from the old and new tables themselves
    awk -F'\t' -v pci="${pci}" 'FILENAME=="ukeys.tmp"{u[$1]=1; next} FILENAME=="o.all"{k=$1"|"$4; if(k in u) po[k]=$pci+0; next} {k=$1"|"$4; if(k in u) pn[k]=$pci+0}
         END{ n=0; for(k in u){ n++; if((k in po) && (k in pn)){ if(pn[k]>po[k]) up++; else if(pn[k]<po[k]) down++ } }
              print "  unobservable burden cells (Number_rare > 0) accepted at region level: " n " (p rose in " up+0 ", fell in " down+0 ")" }' ukeys.tmp o.all n.all >&2
    rm -f ukeys.tmp
  fi
  if [[ ${nbad} -gt 0 ]]; then echo "ERROR ${sfx}: ${nbad} changed row(s) are neither a recorded failure nor a MAC <= ${ERMAC} cell/marker nor an exact-by-weight cell nor a pinned-sample test nor a statistic-identical pseudo-marker row -- splice refused" >&2; exit 1; fi
  echo "${sfx}: $(wc -l < changed.rows) changed row(s), all explained (failed saddlepoint, MAC <= ${ERMAC}, exact test by weight, pinned samples left out, or statistic identity)"
  # CONTROL 2: every recorded failure must have changed. The recorded failure
  # is the PSEUDO-MARKER's own single test (the scan reads the singleAssoc),
  # so it is checked exactly in the singleAssoc; in the results file the
  # burden cell of a set with rare, non-collapsed markers (Number_rare > 0)
  # is a different test that need not have failed, so there it is a note.
  if [[ ${sfx} == ".txt.gz" && -s failcells.txt ]]; then
    nmiss=$(awk -F'\t' 'NR==FNR{ch[$1"|"$4]=1; next} !($1 in ch){m++; if(m<=5) print "  failed pseudo-marker whose burden cell is unchanged: " $1 > "/dev/stderr"} END{print m+0}' changed.rows failcells.txt)
    [[ ${nmiss} -eq 0 ]] || echo "NOTE ${sfx}: ${nmiss} recorded pseudo-marker failure(s) sit in a burden cell that did not change (a set with Number_rare > 0 is a different test; the pseudo-marker row itself is checked in the singleAssoc)" >&2
  fi
  if (( single )) && [[ -s failcells.txt || -s failmarkers.txt ]]; then
    nmiss=$(awk -F'\t' 'NR==FNR{ if('"${isreal}"'){ ch[$3]=1 } else { n=split($3,a,":"); ch[a[1]"|"a[n]]=1 } next} !($1 in ch){m++; if(m<=5) print "  recorded failure unchanged in the singleAssoc: " $1 > "/dev/stderr"} END{print m+0}' changed.rows failcells.txt failmarkers.txt)
    [[ ${nmiss} -eq 0 ]] || echo "WARNING ${sfx}: ${nmiss} recorded failure(s) came back identical in the singleAssoc" >&2
  fi
  # ROW ACCOUNTING: a region-path row whose p is NA is not written (LEDGER
  # #127's convention; since #164 a saddlepoint failure above 13 carriers
  # yields such a row), so a rerun region can come back with FEWER rows than
  # it had. Print every region whose count changed, so the drop is visible.
  # Row counts per region, old against new, in ONE pass per file. This used to
  # be two full passes per REGION, which a filtered rerun of tens of regions
  # never noticed and a complete rerun (REGIONS=ALL, ~18,000 regions against
  # million-row files) would have spent a day on (2026-09-29).
  local why
  if [[ ${sfx} == ".pooledTests.txt.gz" ]]; then why="a region with nothing to disclose -- plain Cauchy, no refused cell -- writes no pooled row (LEDGER #127)"
  else why="NA rows are not written on the region path; a WHOLE region missing is a lost chromosome, read chrStatus (LEDGER #170)"; fi
  awk -F'\t' -v sfx="${sfx}" -v why="${why}" '
    FILENAME=="regions.txt"{keep[$1]=1; next}
    FNR==1{next}
    { reg=('"${rule}"'); if(!(reg in keep)) next; if(FILENAME=="old.tsv") o[reg]++; else n[reg]++ }
    END{ for(r in keep) if((o[r]+0) != (n[r]+0)) print "ROWS " sfx ": " r " old " o[r]+0 " new " n[r]+0 " (" why ")" }' regions.txt old.tsv new.tsv | sort
  if (( single )); then
    awk -F'\t' 'NR==FNR{nm[$1]=$2; next} FNR>1 && '"${isreal}"' && ($3 in nm){ seen[$3]++; if(seen[$3] <= nm[$3]) rep++; else kept++ } END{print rep+0"\t"kept+0}' newmarkers.txt old.tsv > rm.tmp
    echo "${sfx}: real-marker rows: $(wc -l < newmarkers.txt) rerun marker(s), $(cut -f1 rm.tmp) old row(s) replaced by $(awk -F'\t' 'FNR>1 && '"${isreal}"'' new.tsv | wc -l) new, $(cut -f2 rm.tmp) old cop(y/ies) kept for regions not rerun"
  fi
  # the splice: old rows not replaced, then new rows, in old-file order of first appearance
  { head -1 old.tsv
    awk -F'\t' -v rule="${rule}" -v single="${single}" 'FILENAME=="regions.txt"{keep[$1]=1; next} FILENAME=="newmarkers.txt"{nm[$1]=$2; next} FNR==1{next} { if(single && '"${isreal}"'){ if(!($3 in nm)) print; else { seen[$3]++; if(seen[$3] > nm[$3]) print }; next } reg=('"${rule}"'); if(!(reg in keep)) print }' regions.txt newmarkers.txt old.tsv
    awk -F'\t' 'FNR>1' new.tsv
  } > "${out}"
  gzip -f "${out}"
  echo "${sfx}: old $(($(wc -l < old.tsv)-1)) rows, new $(($(wc -l < new.tsv)-1)) rows, spliced $(( $(gzip -dc "${out}.gz" | wc -l) - 1 )) rows"
}
splice_one ".txt.gz"             '$1'
if exists "${OLD_PREFIX}.singleAssoc.txt.gz"; then
  splice_one ".singleAssoc.txt.gz" 'substr($3, 1, index($3, ":") ? index($3, ":")-1 : length($3))'
fi
if exists "${OLD_PREFIX}.pooledTests.txt.gz"; then
  splice_one ".pooledTests.txt.gz" '$1'
fi
if exists "${NEW_PREFIX}.exactByWeight.txt.gz"; then   # LEDGER #172: the routed cells of the spliced regions travel with the file
  # one pass that reads the whole stream: `| head -1` on a list of 1e5 rows closes
  # the pipe under the writer, and pipefail makes that SIGPIPE fatal (exit 141)
  cat_any "${NEW_PREFIX}.exactByWeight.txt.gz" | awk -F'\t' 'NR==FNR{keep[$1]=1; next} FNR==1{print; next} ($2 in keep)' regions.txt - > "${OUT_DIR}/${pre}.exactByWeight.txt"
  gzip -f "${OUT_DIR}/${pre}.exactByWeight.txt"
  echo ".exactByWeight.txt.gz: $(( $(gzip -dc "${OUT_DIR}/${pre}.exactByWeight.txt.gz" | wc -l) - 1 )) routed cell(s) in the spliced regions, carried over"
fi
if exists "${NEW_PREFIX}.stretchGate.txt.gz"; then   # LEDGER #176: the gated cells of the spliced regions travel with the file
  cat_any "${NEW_PREFIX}.stretchGate.txt.gz" | awk -F'\t' 'NR==FNR{keep[$1]=1; next} FNR==1{print; next} ($3 in keep)' regions.txt - > "${OUT_DIR}/${pre}.stretchGate.txt"
  gzip -f "${OUT_DIR}/${pre}.stretchGate.txt"
  echo ".stretchGate.txt.gz: $(gzip -dc "${OUT_DIR}/${pre}.stretchGate.txt.gz" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; c=("reported" in h)?h["reported"]:11; next} (NF>=c && $c=="exact"){n++} END{print n+0}') gated cell(s) in the spliced regions, carried over"
fi
if exists "${NEW_PREFIX}.spaPinned.txt.gz"; then   # LEDGER #174: the guarded tests of the spliced regions travel with the file
  cat_any "${NEW_PREFIX}.spaPinned.txt.gz" | awk -F'\t' 'NR==FNR{keep[$1]=1; next} FNR==1{print; next} ($3 in keep)' regions.txt - > "${OUT_DIR}/${pre}.spaPinned.txt"
  gzip -f "${OUT_DIR}/${pre}.spaPinned.txt"
  echo ".spaPinned.txt.gz: $(gzip -dc "${OUT_DIR}/${pre}.spaPinned.txt.gz" | awk 'END{print NR-1}') guarded test(s) in the spliced regions, carried over"
fi
{ echo "spliced $(date -u +%Y-%m-%dT%H:%M:%SZ)"; echo "old: ${OLD_PREFIX}"; echo "new: ${NEW_PREFIX}"
  echo "regions replaced: $(wc -l < regions.txt)"; echo "regions with failed saddlepoints: $(wc -l < failregions.txt)"
  echo "rows accepted on statistic identity alone (CONTROL 1c, the weaker tier): $( [[ -f "${OUT_DIR}/${pre}.spliceStatId.txt" ]] && awk 'END{print NR}' "${OUT_DIR}/${pre}.spliceStatId.txt" || echo 0 )"
  echo "--- regions:"; cat regions.txt; } > "${OUT_DIR}/${pre}.splice.txt"
echo done
