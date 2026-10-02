import numpy as np, math
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spa_exact_ref import reference, pmf_grid, tail

rng = np.random.default_rng(1)
ok = []

# ---- CONTROL 1 (positive): equal weights -> S is Binomial, closed form -----
n, p = 40, 0.1122
g = np.r_[np.ones(n), np.zeros(5)]; mu = np.full(n + 5, p); nz = np.r_[np.ones(n), np.zeros(5)]
for k in (20, 30, 38):
    exact = sum(math.comb(n, j) * p**j * (1-p)**(n-j) for j in range(k, n+1))
    r = reference(g, mu, nz, float(k) - 1e-9, -1e30, verbose=False)
    lo, hi = r["P(S>=q)"]
    # the bracket is degenerate here (weights land exactly on the lattice), so
    # compare with a relative tolerance rather than demanding containment to
    # the last bit -- the closed form and the DP accumulate rounding differently
    good = abs(hi - exact) <= 1e-9 * exact and lo <= hi
    ok.append(good)
    print(f"{'ok  ' if good else 'FAIL'} binomial tail k>={k:3d}: exact {exact:.6e}  bracket [{lo:.6e}, {hi:.6e}]")

# ---- CONTROL 2 (negative): q at the mean -> p ~ 0.5 ------------------------
g = rng.random(60); mu = rng.uniform(.02, .4, 60)
nz = np.ones(60)
q = float(np.sum(g * mu))
r = reference(g, mu, nz, q, -1e30, verbose=False)
lo, hi = r["P(S>=q)"]
good = 0.4 < lo < 0.6
ok.append(good); print(f"{'ok  ' if good else 'FAIL'} null: q at the mean -> P(S>=q) = [{lo:.4f}, {hi:.4f}] (want ~0.5)")

# ---- CONTROL 3: agrees with Monte Carlo WHERE MC HAS POWER -----------------
# one atom + many tiny, but q chosen so p ~ 1e-3
g = np.r_[1.0, np.full(81, 0.0027)]; mu = rng.uniform(.012, .49, 82); nzf = np.ones(82)
gNA = rng.normal(0, .002, 2000); muNA = rng.uniform(.02, .4, 2000)
G = np.r_[g, gNA]; MU = np.r_[mu, muNA]; NZ = np.r_[nzf, np.zeros(2000)]
rm = float(np.sum(gNA*muNA)); rs = float(np.sqrt(np.sum(muNA*(1-muNA)*gNA**2)))
Y = rng.binomial(1, mu[:, None], size=(82, 4_000_000))
S = (Y * g[:, None]).sum(0) + rng.normal(rm, rs, 4_000_000)
for want in (1e-2, 1e-3):
    qq = float(np.quantile(S, 1 - want))
    mc = float((S >= qq).mean()); se = math.sqrt(mc * (1 - mc) / len(S))
    lo, hi = reference(G, MU, NZ, qq, -1e30, verbose=False)["P(S>=q)"]
    good = lo - 4*se <= mc <= hi + 4*se
    ok.append(good)
    print(f"{'ok  ' if good else 'FAIL'} vs MC at p~{want:g}: MC {mc:.4e} (+-{se:.1e})  conv [{lo:.4e}, {hi:.4e}]")

# ---- CONTROL 4: the atom-dominated regime ---------------------------------
qq = 0.999 * 1.0 + rm          # just under the atom alone
lo, hi = reference(G, MU, NZ, qq, -1e30, verbose=False)["P(S>=q)"]
print(f"\natom-dominated probe: q just under the atom's own weight -> P = [{lo:.4g}, {hi:.4g}]"
      f"   (mu of the atom carrier = {mu[0]:.4f})")


# ===================== PHASE 2: does Fix 1 leak? ==========================
# Fix 1 = enumerate the heaviest 13 carriers exactly, fold the rest into the
# Gaussian remainder. Oracle = exact convolution over all of them. Both are
# given the SAME q, so this isolates the approximation, not the q defect.
print("\n--- phase 2: Fix 1 (top 13 exact, rest Gaussian) against the oracle ---")
MU, RS, RM = 0.01, 0.001, -0.01        # round SYNTHETIC values in the regime of interest

def _t(w, mu, q, rm, rs):
    step = max(w.sum() / 2.0e5, 1e-13); v = []
    for m in ("down", "up"):
        lo, pmf = pmf_grid(w, mu, step, m); v.append(tail(lo, pmf, step, q, rm, rs, +1))
    return (min(v) + max(v)) / 2

def fix1_vs_oracle(w, q, topm=13):
    mu = np.full(len(w), MU)
    o = _t(w, mu, q, RM, RS)
    i = np.argsort(w)[::-1]; k, d = i[:topm], i[topm:]
    rm2 = RM + float(np.sum(w[d] * mu[d]))
    rs2 = float(np.sqrt(RS**2 + np.sum(w[d]**2 * mu[d] * (1 - mu[d]))))
    f = _t(w[k], mu[k], q, rm2, rs2)
    neff = w.sum() ** 2 / (w ** 2).sum()
    return neff, o, f

# (a) a near-one-atom shape: one heavy weight, many negligible ones
w = np.r_[1.0, np.full(81, 0.003 / 81)]
for label, q in (("one-atom shape, correct q", 0.99),
                 ("one-atom shape, inflated q (the bug)", 0.99 / math.sqrt(0.98))):
    ne, o, f = fix1_vs_oracle(w, q)
    good = 0.95 < f / o < 1.05
    ok.append(good)
    print(f"{'ok  ' if good else 'FAIL'} {label:34s} n_eff {ne:6.3f}  "
          f"oracle {o:.4e}  fix1 {f:.4e}  ratio {f/o:6.3f}")

# (b) adversarial: one atom plus k EQUAL secondary carriers. Equal weights are a
#     plausible adversary, NOT a proven worst case -- real cells vary in both
#     weight and mu. Only shapes the routing would admit (n_eff <= 4) count.
print("\n     adversarial sweep (admitted shapes only; worst case NOT proven):")
worst_dec, worst_any = 1.0, 1.0
for k in (2, 5, 13, 14, 20, 30, 60):
    for T in (0.05, 0.25, 1.0):
        ne, o, f = fix1_vs_oracle(np.r_[1.0, np.full(k, T / k)], 1.0)
        if ne > 4.0 or o <= 0:
            continue
        r = f / o
        worst_any = max(worst_any, r, 1 / r)
        if o >= 1e-3:
            worst_dec = max(worst_dec, r, 1 / r)
print(f"     worst error anywhere: {worst_any:6.1f}x   "
      f"worst at p >= 1e-3 (decision-relevant): {worst_dec:5.2f}x")
print("\nall controls passed:", all(ok))
