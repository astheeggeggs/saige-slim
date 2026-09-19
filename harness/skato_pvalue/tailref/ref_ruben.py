"""Cancellation-free reference for P(sum_j lambda_j chi2_1 > q), small q (4-10 terms).

Ruben/Farebrother expands  F(x) = sum_k a_k P(chi2_{n+2k} <= x/beta).
Since sum_k a_k = 1, the SURVIVAL is  S(x) = sum_k a_k P(chi2_{n+2k} > x/beta).
Choosing beta <= min(lambda) makes every gamma_j = 1 - beta/lambda_j lie in [0,1),
hence every g_r >= 0 and every a_k >= 0: the survival is a sum of NONNEGATIVE
terms and there is no 1-F cancellation.  That is the whole point -- Davies,
Farebrother and Imhof all reach the far tail by subtracting from 1, which is why
they fail there, and it is why this reference can be trusted where they cannot.

Truncating at K omits mass R_K = 1 - sum_{k<=K} a_k, and since every chi2 factor
is <= 1 that is a strict upper bound on the omitted contribution.  A value is
returned as ok only when R_K < rtol * S.  The reference REFUSES rather than
guesses.
"""
import numpy as np
from scipy.stats import chi2


def ruben_survival(lam, q, Kmax=20000, rtol=1e-3):
    """Adaptive: grow the series until the omitted coefficient mass R_K is
    provably below rtol * (partial survival).  R_K = 1 - sum_{k<=K} a_k bounds
    the omitted contribution because every chi2 factor is <= 1."""
    lam = np.asarray(lam, dtype=float)
    n = lam.size
    beta = lam.min() * (1 - 1e-13)
    gam = 1.0 - beta / lam
    x = q / beta
    a = np.empty(Kmax + 1)
    g = np.empty(Kmax + 1)
    a[0] = np.exp(0.5 * np.log(beta / lam).sum())
    S = a[0] * chi2.sf(x, n)
    asum = a[0]
    pw = np.ones_like(gam)
    K = 0
    CHUNK = 64
    for k in range(1, Kmax + 1):
        pw *= gam
        g[k] = pw.sum()
        a[k] = float(np.dot(g[1:k + 1], a[k - 1::-1])) / (2.0 * k)
        S += a[k] * chi2.sf(x, n + 2 * k)
        asum += a[k]
        K = k
        if k % CHUNK == 0:
            R = 1.0 - asum
            if R <= 0 or (S > 0 and R < rtol * S):
                break
    R = max(1.0 - asum, 0.0)
    return S, R, (S > 0 and R < rtol * S), K


# ---- exact closed forms, the positive controls -------------------------------
def exact_equal(lam_value, n, q):
    """lambda_j == c for all j  =>  Q = c * chi2_n, exact."""
    return float(chi2.sf(q / lam_value, n))


def exact_two_pairs(a, b, q):
    """lambda = (a,a,b,b) => Q = 2a E1 + 2b E2, E~Exp(1), a != b.
    S(q) = (alpha e^{-beta q} - beta e^{-alpha q}) / (alpha - beta),
    alpha = 1/(2a), beta = 1/(2b)."""
    al, be = 1.0 / (2 * a), 1.0 / (2 * b)
    return float((al * np.exp(-be * q) - be * np.exp(-al * q)) / (al - be))
