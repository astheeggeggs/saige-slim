#!/usr/bin/env python3
"""Line-by-line comparison of two analyse.R outputs (results-<date>.txt).

usage: compare_results.py OLD NEW
Keys each CONTROL 1 / 2 / 3 and COLLAPSE line by its labels and prints every
line whose numbers differ, grouped by method class, plus a per-class summary
(max |delta obs| at each alpha, max |delta detections|, max |delta rho|).
Written for the 2026-09-12 re-run at HEAD against the 2026-09-10 sweep
(docs/state/prediction-flexrv-rerun-2026-09-12.md: skato/grid must not move
at all; unw/beta by <= 2 genes; flex* by <= 5 genes / 3 detections)."""
import re, sys, collections

def parse(fn):
    sec = None; out = {}
    for line in open(fn):
        if line.startswith("CONTROL 1"): sec = "C1"; continue
        if line.startswith("CONTROL 2"): sec = "C2"; continue
        if line.startswith("CONTROL 3"): sec = "C3"; continue
        if line.startswith("COLLAPSE"): sec = "COL"; continue
        if not line.startswith("  ") or sec is None: continue
        t = line.split()
        if sec == "C1" and "alpha=" in line:
            i = [k for k, x in enumerate(t) if x.startswith("alpha=")][0]
            key = ("C1",) + tuple(t[:i]) + (t[i],)
            m = re.search(r"obs=([0-9.]+)\s+\[([0-9.]+),\s*([0-9.]+)\]\s+n=(\d+)", line)
            out[key] = {"obs": float(m.group(1)), "lo": float(m.group(2)), "hi": float(m.group(3)),
                        "n": int(m.group(4)), "outside": "OUTSIDE" in line, "method": t[i - 2]}
        elif sec == "C2" and "rho=" in line:
            key = ("C2",) + tuple(t[:4])
            m = re.search(r"rho=\s*(-?[0-9.]+)\s+detected\s+(\d+)/(\d+)\s+\(directional\s+(\d+)/(\d+)\)\s+oracle\s+(\S+)", line)
            out[key] = {"rho": float(m.group(1)), "det": int(m.group(2)), "dir": int(m.group(4)),
                        "oracle": m.group(6), "method": t[3], "bylink": line.split("by link")[1].strip()}
        elif sec == "C3" and "flex192:" in line:
            key = ("C3",) + tuple(t[:2])
            nums = [int(x) for x in re.findall(r"(?:informative|realistic|uninformative|unw)\s+(\d+)", line)]
            out[key] = {"counts": nums, "flags": re.findall(r"(monotone=\w+|uninf<=unw=\w+)", line), "method": "flex192"}
        elif sec == "COL" and "median=" in line:
            key = ("COL",) + tuple(t[:2])
            m = re.search(r"median=([0-9.]+)\s+max=([0-9.]+).*either:\s*(\d+)", line)
            out[key] = {"median": float(m.group(1)), "max": float(m.group(2)), "n1e3": int(m.group(3)), "method": "flex192"}
    return out

def cls(method):
    if method in ("skato", "grid"): return "A: skato/grid (must not move)"
    if method in ("unw", "beta"): return "B: unw/beta (<= 2 genes)"
    return "C: flex* (may move)"

old, new = parse(sys.argv[1]), parse(sys.argv[2])
print(f"old lines {len(old)}  new lines {len(new)}  only-old {len(set(old)-set(new))}  only-new {len(set(new)-set(old))}")
summ = collections.defaultdict(lambda: {"lines": 0, "changed": 0, "dobs": 0.0, "ddet": 0, "drho": 0.0, "dgenes": 0.0})
for k in sorted(set(old) & set(new)):
    a, b = old[k], new[k]; c = cls(a["method"]); s = summ[c]; s["lines"] += 1
    if a == b: continue
    s["changed"] += 1
    if k[0] == "C1":
        d = b["obs"] - a["obs"]; s["dobs"] = max(s["dobs"], abs(d)); s["dgenes"] = max(s["dgenes"], abs(d) * a["n"])
        print(f"{c[:1]} {' '.join(k[1:])}: obs {a['obs']:.4f} -> {b['obs']:.4f} ({d * a['n']:+.1f} genes){' OUTSIDE' if b['outside'] else ''}")
    elif k[0] == "C2":
        s["ddet"] = max(s["ddet"], abs(b["det"] - a["det"])); s["drho"] = max(s["drho"], abs(b["rho"] - a["rho"]))
        print(f"{c[:1]} {' '.join(k[1:])}: rho {a['rho']:+.2f} -> {b['rho']:+.2f}, detected {a['det']} -> {b['det']}, oracle {a['oracle']} -> {b['oracle']}"
              + ("" if a["bylink"] == b["bylink"] else f", by link {a['bylink']} -> {b['bylink']}"))
    elif k[0] == "C3":
        print(f"{c[:1]} {' '.join(k[1:])}: counts {a['counts']} -> {b['counts']}, flags {a['flags']} -> {b['flags']}")
    else:
        print(f"{c[:1]} {' '.join(k[1:])}: median {a['median']} -> {b['median']}, max {a['max']} -> {b['max']}, n(p<1e-3) {a['n1e3']} -> {b['n1e3']}")
print("\nsummary by class:")
for c in sorted(summ):
    s = summ[c]
    print(f"  {c}: {s['changed']} of {s['lines']} lines changed; max |d obs| {s['dobs']:.4f} ({s['dgenes']:.1f} genes); max |d detected| {s['ddet']}; max |d rho| {s['drho']:.2f}")
