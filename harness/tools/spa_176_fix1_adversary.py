import sys, os, math, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spa_exact_ref import pmf_grid, tail, reference
rng = np.random.default_rng(7)

# ---------- A. ORACLE vs brute-force enumeration (bracket containment) ----------
def brute(w, mu, q, rm, rs, side):
    n = len(w); tot = 0.0
    for bits in itertools.product((0, 1), repeat=n):
        b = np.array(bits); pr = np.prod(np.where(b, mu, 1 - mu)); s = float(w @ b)
        if rs <= 0:
            tot += pr * ((s + rm >= q) if side > 0 else (s + rm <= q))
        else:
            z = (q - s - rm) / rs
            tot += pr * 0.5 * math.erfc((z if side > 0 else -z) / math.sqrt(2))
    return tot
viol = 0; ntest = 0; widest = 0.0
for trial in range(120):
    n = int(rng.integers(3, 13))
    kind = trial % 4
    if kind == 0: w = rng.random(n)
    elif kind == 1: w = rng.normal(0, 1, n)                 # signed (gtilde can be negative)
    elif kind == 2: w = np.r_[1.0, rng.random(n - 1) * 1e-3] # atom + tiny
    else: w = np.r_[1.0, rng.random(n - 1) * 1e-7]           # below the default step
    mu = 10 ** rng.uniform(-4, -0.3, n)
    rs = [0.0, 1e-3, 0.05][trial % 3]; rm = rng.normal(0, .01)
    full = np.r_[w, 0.0]; MU = np.r_[mu, 0.5]; NZ = np.r_[np.ones(n), 0]
    # reference() derives rm/rs from non-carriers, so call the pieces directly
    span = np.abs(w).sum(); step = max(span / 2e5, 1e-12)
    for q in list(rng.uniform(w[w<0].sum()-.01, w[w>0].sum()+.01, 3)) + [float(w.max())+rm, float(w @ (mu>0.2))+rm]:
        for side in (1, -1):
            ex = brute(w, mu, q, rm, rs, side)
            v = [tail(*pmf_grid(w, mu, step, m), step, q, rm, rs, side) for m in ("down", "up")]
            lo, hi = min(v), max(v); ntest += 1
            if not (lo <= ex * (1 + 1e-9) + 1e-300 and ex <= hi * (1 + 1e-9) + 1e-300):
                viol += 1
                if viol <= 8: print(f"  VIOLATION kind {kind} rs {rs} side {side:+d} q {q:.6g}: exact {ex:.6e} bracket [{lo:.6e}, {hi:.6e}]")
            if ex > 1e-300 and lo > 0: widest = max(widest, hi / lo)
            elif hi > 0 and lo == 0: widest = float("inf")
print(f"A. bracket containment vs 2^n enumeration: {viol} violations in {ntest}; widest hi/lo = {widest}")

# ---------- B. bracket WIDTH in phase 2's own cases (it reports midpoints) ----------
MU0, RS, RM = 0.01, 0.001, -0.01      # round SYNTHETIC values
def br(w, mu, q, rm, rs):
    step = max(w.sum() / 2e5, 1e-13)
    v = [tail(*pmf_grid(w, mu, step, m), step, q, rm, rs, +1) for m in ("down", "up")]
    return min(v), max(v)
w = np.r_[1.0, np.full(81, 0.003 / 81)]; mu = np.full(82, MU0)
for lab, q in (("correct q", 0.99), ("inflated q", 0.99 / math.sqrt(0.98))):
    lo, hi = br(w, mu, q, RM, RS); print(f"B. oracle bracket, one-atom shape, {lab}: [{lo:.4e}, {hi:.4e}]  hi/lo {hi/lo:.3f}")

# ---------- C. Fix 1 adversary ----------
def fix1(w, mu, q, rm, rs, topm=13):
    i = np.argsort(w)[::-1]; k, d = i[:topm], i[topm:]
    rm2 = rm + float(np.sum(w[d] * mu[d])); rs2 = math.sqrt(rs**2 + float(np.sum(w[d]**2 * mu[d] * (1 - mu[d]))))
    return br(w[k], mu[k], q, rm2, rs2)
def case(lab, w, mu, q, rm=RM, rs=RS):
    neff = w.sum()**2 / (w**2).sum()
    olo, ohi = br(w, mu, q, rm, rs); flo, fhi = fix1(w, mu, q, rm, rs)
    mean = float(w @ mu) + rm; sd = math.sqrt(float(np.sum(w*w*mu*(1-mu))) + rs*rs)
    o = math.sqrt(olo*ohi); f = math.sqrt(max(flo,1e-320)*max(fhi,1e-320))
    return dict(lab=lab, neff=neff, z=(q-mean)/sd, o=o, f=f, r=f/o, ow=ohi/olo)
rows = []
# family 1 (theirs, for baseline) : equal secondaries, single mu, q = 1.0
# family 2: atom + k equal secondaries, LOW mu on the secondaries, q at ATTAINABLE observed points:
#           atom fired (+ j secondaries fired), or atom NOT fired and j secondaries fired
for mua in (0.01, 0.1, 0.3):
  for mus in (0.01, 0.003, 0.001):
    for k in (20, 40, 80, 200):
      for T in (0.25, 0.5, 1.0):
        a = T / k
        w = np.r_[1.0, np.full(k, a)]; mu = np.r_[mua, np.full(k, mus)]
        if w.sum()**2/(w**2).sum() > 4.0: continue
        base = RM   # remainder at its mean
        for j in (1, 2, 3, 4, 6):
            rows.append(case(f"atom FIRED +{j} sec | mua {mua} mus {mus} k {k} T {T}", w, mu, 1.0 + j*a + base))
            rows.append(case(f"atom not fired, {j} sec | mua {mua} mus {mus} k {k} T {T}", w, mu, j*a + base))
# family 3: geometric secondaries with heterogeneous mu
for mus in (0.003, 0.001):
  for rho in (0.9, 0.95, 0.98):
    k = 80; a = 0.08 * rho ** np.arange(k); w = np.r_[1.0, a]; mu = np.r_[0.05, np.full(k, mus)]
    if w.sum()**2/(w**2).sum() > 4.0: continue
    for j in (1, 2, 3):
        rows.append(case(f"geom rho {rho} mus {mus}: atom not fired, carriers 14..{13+j} fired", w, mu, a[13:13+j].sum() + RM))
        rows.append(case(f"geom rho {rho} mus {mus}: atom FIRED + carriers 14..{13+j}", w, mu, 1.0 + a[13:13+j].sum() + RM))
rows = [r for r in rows if r["z"] > 2 and r["o"] > 0]
print(f"\nC. {len(rows)} admitted (n_eff<=4), z>2 probes at ATTAINABLE q. Worst fix1/oracle by true-p decade (ratio<1 = fix1 ANTI-conservative):")
for lo_, hi_ in ((1e-2, 1), (1e-3, 1e-2), (1e-4, 1e-3), (1e-6, 1e-4), (1e-8, 1e-6), (1e-12, 1e-8)):
    sel = [r for r in rows if lo_ <= r["o"] < hi_]
    if not sel: continue
    mn = min(sel, key=lambda r: r["r"]); mx = max(sel, key=lambda r: r["r"])
    print(f"  p in [{lo_:g},{hi_:g}) n={len(sel):4d}  min ratio {mn['r']:.3g} (o {mn['o']:.2e} f {mn['f']:.2e} neff {mn['neff']:.2f} z {mn['z']:.1f}; {mn['lab']})")
    print(f"  {'':27s}  max ratio {mx['r']:.3g} (o {mx['o']:.2e} f {mx['f']:.2e} neff {mx['neff']:.2f}; {mx['lab']})")
    print(f"  {'':27s}  share off by >2x: {np.mean([abs(math.log(r['r']))>math.log(2) for r in sel]):.2f}; oracle bracket hi/lo max {max(r['ow'] for r in sel):.2f}")
