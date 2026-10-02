import sys, os, math, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spa_exact_ref import pmf_grid, tail
VR = 0.98; MU_ATOM = 0.01; M1 = 0.004   # round SYNTHETIC values, no participant's

def Phibar(z): return 0.5 * math.erfc(z / math.sqrt(2))

def spa_tail(w, mu, rm, rs, q):
    """Lugannani-Rice in SAIGE's Barndorff-Nielsen form: P(S >= q), q above the mean."""
    K  = lambda t: float(np.sum(np.logaddexp(np.log1p(-mu), np.log(mu) + w * t))) + rm * t + 0.5 * rs * rs * t * t
    K1 = lambda t: float(np.sum(w * mu / (mu + (1 - mu) * np.exp(-w * t)))) + rm + rs * rs * t
    def K2(t):
        e = np.exp(-w * t); return float(np.sum(w * w * mu * (1 - mu) * e / (mu + (1 - mu) * e) ** 2)) + rs * rs
    lo, hi = 0.0, 1.0
    while K1(hi) < q:
        hi *= 2
        if hi > 1e9: return float("nan")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if K1(mid) < q: lo = mid
        else: hi = mid
    t = 0.5 * (lo + hi)
    ww = math.sqrt(max(2 * (t * q - K(t)), 0)); v = t * math.sqrt(K2(t))
    if ww <= 0 or v <= 0: return float('nan')
    return Phibar(ww + math.log(v / ww) / ww)

def exact_tail(w, mu, rm, rs, q):
    step = max(np.abs(w).sum() / 4e5, 1e-13)
    v = [tail(*pmf_grid(w, mu, step, m), step, q, rm, rs, +1) for m in ("down", "up")]
    return math.sqrt(min(v) * max(v)), max(v) / min(v)

def cell(T, nasigma, shape, seed, n=82):
    rng = np.random.default_rng(seed)
    rest = T
    raw = {"equal": np.ones(n - 1), "geom": 0.93 ** np.arange(n - 1),
           "lognorm": rng.lognormal(0, 1.0, n - 1)}[shape]
    sec = np.minimum(raw / raw.sum() * rest, 0.999)  # no secondary may outweigh the atom
    w = np.r_[1.0, sec]
    mu = np.r_[MU_ATOM, 10 ** rng.uniform(math.log10(5e-3), math.log10(0.487), n - 1)]
    rs = math.sqrt(nasigma); rm = M1 - float(w @ mu)
    sd_rest = math.sqrt(float(np.sum(sec ** 2 * mu[1:] * (1 - mu[1:]))) + rs * rs)
    return w, mu, rm, rs, sd_rest


from statistics import NormalDist
Q_STR = 1.0; Q_OBS = M1 + (Q_STR - M1) * math.sqrt(VR); DELTA = Q_STR - Q_OBS
print(f"stretch: q_obs {Q_OBS:.5f} -> q {Q_STR:.5f} (delta {DELTA:.5f})")
print("T = total weight of the 81 secondary carriers (atom = 1)\n")
print("   T     shape    n_eff  d/sd_rest | SPA(q)     exact(q)   SPA/exact | truth=exact(q_obs)  reported/truth | VR-on-p-scale/truth")
ctrl = {0.003: "real cells of this kind: saddlepoint collapses by many orders",
        0.22:  "real cells of this kind: saddlepoint several-fold conservative at the same q"}
for T, nas in ((0.003, 1.36e-6), (0.01, 1.4e-6), (0.03, 1.5e-6), (0.1, 1.6e-6), (0.22, 1.7e-6),
               (0.5, 5e-6), (1.0, 1.5e-5), (2.3, 3.2e-5), (5.0, 1.6e-4)):
    for shape in ("equal", "geom", "lognorm"):
        w, mu, rm, rs, sdr = cell(T, nas, shape, seed=11)
        neff = w.sum() ** 2 / (w ** 2).sum()
        sp = spa_tail(w, mu, rm, rs, Q_STR); e, bw = exact_tail(w, mu, rm, rs, Q_STR)
        truth, _ = exact_tail(w, mu, rm, rs, Q_OBS)
        z = NormalDist().inv_cdf(1 - truth); pfix = Phibar(z / math.sqrt(VR))
        print(f"{T:6.3f}  {shape:8s} {neff:5.2f} {DELTA/sdr:9.2f}  | {sp:.3e}  {e:.3e}  {sp/e:8.3g}  | {truth:.3e}          {sp/truth:10.3g}     | {pfix/truth:6.3f}")
    if T in ctrl: print("        control ->", ctrl[T])
