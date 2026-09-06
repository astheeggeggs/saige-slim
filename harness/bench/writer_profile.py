#!/usr/bin/env python3
"""Read macOS `sample` call trees from harness/profile.sh and report, per file,
where the MAIN thread's time inside mainMarkerInCPP went.

    python3 harness/bench/writer_profile.py harness/bench/profiles/writer_*.txt

Per profile: main-thread samples inside the mainMarkerInCPP wrapper, split
into
  writer   writeOutfile_single and everything under it
  scan     the serial marker loop (T = 1) or the main thread's own share of
           the parallelFor work (TBB runs tasks on the calling thread too)
  wait     the rest of the parallelFor subtree: TBB scheduling and waiting
  other    what remains (results allocation, post-passes, progress lines)
and the worker threads' totals with their top symbols. Each sample is ~1 ms
of that thread's wall time, so writer/total on the main thread is w.

A node is counted once, into the category of its NEAREST categorised
ancestor-or-self, so nested matches never double count."""
import re
import sys
from collections import Counter

NODE = re.compile(r'^(?P<prefix>[\s+!:|]*?)(?P<n>\d+) (?P<sym>.*?)(?:  \(in |\s+\[0x|$)')
THREAD = re.compile(r'^\s+(\d+) Thread_\d+')

CATS = [
    ('writer', re.compile(r'writeOutfile_single')),
    ('scan', re.compile(r'testOneMarker|prepareOneMarker|testOneModel|MarkerScanWorker|'
                        r'getMarkerPval|Unified_getOneMarker|getOneMarker|scoreTest|'
                        r'imputeGenoAndFlip|countOneMarker')),
    ('wait', re.compile(r'parallelFor|parallel_for|tbb::|wait_for_all|spin|rml::|'
                        r'__psynch|pthread_cond|semaphore|Ppark')),
]
KERNELS = re.compile(r'fast_logistf|Firth|SPA|saddle|getMarkerPval|scoreTest|getOneMarker|'
                     r'countOneMarker|imputeGenoAndFlip|Binary_|ER_|ComputeExact|dgemv|dgemm|'
                     r'daxpy|ddot|arma::|elem\(|Unified_')


def parse(path):
    threads = []   # (header_n, is_main, [(depth, n, sym)])
    cur = None
    for line in open(path, errors='replace'):
        m = THREAD.match(line)
        if m:
            cur = [int(m.group(1)), 'main-thread' in line, []]
            threads.append(cur)
            continue
        if cur is None:
            continue
        m = NODE.match(line.rstrip('\n'))
        if not m:
            continue
        cur[2].append((len(m.group('prefix')), int(m.group('n')), m.group('sym')))
    return threads


def categorise(nodes, root_re):
    """Sum samples under every node matching root_re, by category."""
    out = Counter()
    root_total = 0
    i = 0
    while i < len(nodes):
        d, n, sym = nodes[i]
        if not root_re.search(sym):
            i += 1
            continue
        root_total += n
        stack = []   # (depth, category)
        j = i + 1
        while j < len(nodes) and nodes[j][0] > d:
            dj, nj, sj = nodes[j]
            while stack and stack[-1][0] >= dj:
                stack.pop()
            cat = next((c for c, rx in CATS if rx.search(sj)), None)
            inherited = stack[-1][1] if stack else None
            if cat is not None and inherited is None:
                out[cat] += nj
                stack.append((dj, cat))
            elif cat == 'scan' and inherited == 'wait':
                # TBB runs tasks on the calling thread too: marker work under
                # the parallelFor frame is the main thread's scan share, not
                # waiting. Move it out of the wait bucket it inherited.
                out['wait'] -= nj
                out['scan'] += nj
                stack.append((dj, 'scan'))
            else:
                stack.append((dj, inherited))
            j += 1
        i = j
    out['other'] = root_total - sum(out.values())
    return root_total, out


def worker_summary(threads):
    total = 0
    tops = Counter()
    for hdr, is_main, nodes in threads:
        if is_main:
            continue
        total += hdr
        seen = set()
        for d, n, sym in nodes:
            key = sym.split('(')[0]
            if KERNELS.search(sym) and key not in seen:
                tops[key] += n
                seen.add(key)   # count each symbol once per thread (its top frame)
    return total, tops


def main(paths):
    for p in paths:
        threads = parse(p)
        main_t = [t for t in threads if t[1]]
        if not main_t:
            print(f'{p}: no main thread found')
            continue
        hdr, _, nodes = main_t[0]
        total, cats = categorise(nodes, re.compile(r'_SAIGE_mainMarkerInCPP|^mainMarkerInCPP\('))
        wtotal, tops = worker_summary(threads)
        print(f'== {p}')
        print(f'   main thread samples: {hdr}; inside mainMarkerInCPP: {total} '
              f'({100.0 * total / max(hdr, 1):.1f}% of main)')
        for c in ('writer', 'scan', 'wait', 'other'):
            print(f'   {c:7s} {cats[c]:8d}  {100.0 * cats[c] / max(total, 1):5.1f}% of mainMarkerInCPP')
        print(f'   w = writer / mainMarkerInCPP = {cats["writer"] / max(total, 1):.4f}')
        print(f'   worker threads: {len(threads) - 1} others, {wtotal} samples')
        for sym, n in tops.most_common(8):
            print(f'      {n:8d}  {sym[:90]}')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
