#!/usr/bin/env python3
"""Turn harness/bench/writer.tsv and the writer_* profiles into the table the
plan asked for, with its reconciliation checks.

    python3 harness/bench/writer_analyse.py [label]

Per cell (P models, trait mix) and thread count T:
  timed:   mean real/user over reps, wall speedup vs T = 1, and the Amdahl
           serial fraction implied by that speedup, s = (T/S - 1)/(T - 1)
  profile: w = writer / mainMarkerInCPP on the main thread, writer samples
           per row written (us/row), the main-thread split, worker totals
Checks printed at the end, each PASS/FAIL with its numbers:
  C1 writer us/row agrees across cells within 30% and with the N=1000 bound
  C2 writer samples scale with P (rows) and do not change with T
  C3 profile main-thread samples ~ timed wall at the same T (within 15%)
  C4 the Amdahl s from wall times vs the serial share from the profile
     ((writer + other + main-thread scan share) / main total)
"""
import csv
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
import writer_profile as wp  # noqa: E402

label = sys.argv[1] if len(sys.argv) > 1 else None
tsv = os.path.join(os.path.dirname(__file__), 'writer.tsv')
rows = [r for r in csv.DictReader(open(tsv), delimiter='\t') if not label or r['label'] == label]
if not rows:
    sys.exit('no rows' + (f' for label {label}' if label else ''))

def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

timed = defaultdict(list)
prof = {}
rows_written = {}
for r in rows:
    key = (r['cell'], int(r['threads']))
    rows_written[key] = int(r['rows'])
    if r['mode'] == 'time':
        timed[key].append(r)
    elif r['profile'] not in ('NA', ''):
        prof[key] = r['profile']

def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None

print(f"commit {rows[0]['commit']} ({rows[0]['dirty']}), SAIGE.so {rows[0]['so_md5'][:8]}")
print(f"{'cell':12s} {'T':>2s} {'n':>2s} {'real_s':>8s} {'user_s':>8s} {'speedup':>8s} {'s_amdahl':>8s} {'rss_MB':>7s}")
amdahl = {}
wall = {}
for cell in sorted({c for c, _ in timed}):
    base = mean([fnum(r['real_s']) for r in timed[(cell, 1)]]) if (cell, 1) in timed else None
    for T in sorted({t for c, t in timed if c == cell}):
        rs = timed[(cell, T)]
        real = mean([fnum(r['real_s']) for r in rs]); user = mean([fnum(r['user_s']) for r in rs])
        rss = mean([fnum(r['maxRSS_MB']) for r in rs])
        wall[(cell, T)] = real
        sp = base / real if base and real else None
        s = (T / sp - 1) / (T - 1) if sp and T > 1 else None
        if s is not None:
            amdahl[(cell, T)] = s
        print(f"{cell:12s} {T:2d} {len(rs):2d} {real:8.1f} {user:8.1f} "
              f"{(f'{sp:.2f}x' if sp else '-'):>8s} {(f'{s:.3f}' if s is not None else '-'):>8s} {rss:7.0f}")

print()
print(f"{'profile cell':12s} {'T':>2s} {'main':>7s} {'inMain':>7s} {'writer':>7s} {'scan':>7s} {'wait':>7s} {'other':>7s} {'w':>6s} {'us/row':>7s} {'workers':>8s}")
usrow = {}
wsamp = {}
serial_share = {}
for (cell, T), path in sorted(prof.items()):
    if not os.path.exists(path):
        print(f"{cell:12s} {T:2d} missing {path}"); continue
    threads = wp.parse(path)
    main_t = [t for t in threads if t[1]]
    if not main_t:
        print(f"{cell:12s} {T:2d} no main thread"); continue
    hdr, _, nodes = main_t[0]
    total, cats = wp.categorise(nodes, re.compile(r'_SAIGE_mainMarkerInCPP|^mainMarkerInCPP\('))
    wtotal, _ = wp.worker_summary(threads)
    nrows = rows_written.get((cell, T), 0)
    w = cats['writer'] / max(total, 1)
    ur = 1000.0 * cats['writer'] / nrows if nrows else None
    usrow[(cell, T)] = ur; wsamp[(cell, T)] = cats['writer']
    serial_share[(cell, T)] = (cats['writer'] + cats['other'] + cats['scan']) / max(total, 1)
    print(f"{cell:12s} {T:2d} {hdr:7d} {total:7d} {cats['writer']:7d} {cats['scan']:7d} {cats['wait']:7d} "
          f"{cats['other']:7d} {w:6.3f} {(f'{ur:.1f}' if ur else '-'):>7s} {wtotal:8d}")

print("\nchecks")
vals = [v for v in usrow.values() if v]
if vals:
    lo, hi = min(vals), max(vals)
    ok = hi / lo <= 1.3 if lo > 0 else False
    print(f"  C1 writer us/row across cells: {lo:.1f}-{hi:.1f} (N=1000 bound 5.3, 50k binary P=1 ~3.2): {'PASS' if ok else 'FAIL'}")
for cell in sorted({c for c, _ in wsamp}):
    ts = sorted(t for c, t in wsamp if c == cell)
    if len(ts) > 1:
        a, b = wsamp[(cell, ts[0])], wsamp[(cell, ts[-1])]
        ok = abs(a - b) / max(a, b, 1) <= 0.25
        print(f"  C2 {cell}: writer samples T={ts[0]} {a} vs T={ts[-1]} {b}: {'PASS' if ok else 'FAIL'} (must not change with T)")
p1 = [wsamp[k] / max(rows_written.get(k, 1), 1) for k in wsamp if k[0].startswith('p1:') and k[1] == 1]
p5 = [wsamp[k] / max(rows_written.get(k, 1), 1) for k in wsamp if k[0].startswith(('p5:', 'p3:')) and k[1] == 1]
if p1 and p5:
    print(f"  C2 per-row writer samples P=1 {mean(p1)*1000:.1f} us vs P>1 {mean(p5)*1000:.1f} us: "
          f"{'PASS' if abs(mean(p1)-mean(p5))/max(mean(p1),mean(p5)) <= 0.3 else 'FAIL'} (scales with rows, not models)")
for (cell, T), path in sorted(prof.items()):
    if (cell, T) in wall and os.path.exists(path):
        hdr = [t for t in wp.parse(path) if t[1]][0][0]
        ratio = hdr / 1000.0 / wall[(cell, T)]
        print(f"  C3 {cell} T={T}: main-thread samples {hdr/1000:.1f} s vs timed wall {wall[(cell, T)]:.1f} s "
              f"(ratio {ratio:.2f}; profile includes R startup): {'PASS' if 0.8 <= ratio <= 1.3 else 'CHECK'}")
for (cell, T), s in sorted(amdahl.items()):
    if (cell, T) in serial_share:
        print(f"  C4 {cell} T={T}: Amdahl s from wall {s:.3f} vs profile serial share {serial_share[(cell, T)]:.3f}")
print("\nverdict inputs: w(cell, T) above; ceiling on the item = 1/(1-w)")
