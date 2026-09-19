# The outer integral validated, and the escalation lever costed (2026-09-19)

Two loose ends closed. Neither changes any code.

## 1. The outer rho integral does NOT amplify the inner bias

`handoff-2026-09-19-saddlepoint-tail-validation.md` validated the INNER
survival function (`skatoSaddlepointSurvival`: 3.6-12% conservative, flat in
the number of terms) and explicitly left open how that bias composes through
the direct piecewise outer integral. **Measured: it does not compose. The
outer integral contributes <=0.5% at p ~ 1e-4 and <=0.01% below 1e-8.**

Method (`harness/skato_pvalue/tailref/outer.R`): run the shipped
`skatoSaddlepointTail` against a copy that is identical in every respect --
same rho envelope, same kink-finding, same piecewise quadrature, same
post-rules -- with ONLY the outer integrand's survival call swapped for the
cancellation-free Ruben reference. **The per-rho layer is identical in both
arms**, so `pmin`, `pmin.q` and the envelope are shared and any difference is
the outer integrand's survival propagating. That is exactly the open question,
and isolating it this way is also ~9x cheaper than swapping both layers
(Ruben's coefficient recursion is per-lambda; the outer integrand hits ONE
fixed lambda at every quadrature node, so it is cached once per cell).

Real multi-rho cells from the finished-build guard sweep, at their own depth:

| cells | q | ratio shipped/reference |
|---|---|---|
| 6 unique | 2-39 | 0.99187 .. 1.00448, median 1.00176 (max abs log10 0.0035) |

Deeper than the null reaches, by scaling `Score` (which scales Q, i.e. exactly
what stronger signal looks like, leaving Phi's structure untouched):

| p | ratio |
|---|---|
| 1.5e-04 | 1.00448 |
| 1.1e-14 | 1.00008 |
| 7.4e-54 | 1.00000 |
| 1.6e-211 | 1.00000 |

**The mechanism, checked rather than assumed.** The reported p is
`pchisq(xstar, df = 1, lower.tail = FALSE) + area`, and the first term does not
involve the survival function at all. As the tail deepens it takes over:
92.0% of the p at 1.5e-4, 97.8% at 6.0e-7, 99.3% at 6.9e-10, 99.8% at 1.1e-14.
So the outer integrand's error is diluted precisely where the tail method fires.

**The vacuity check that had to be run first.** `skatoSaddlepointTail` caps at
`pmin * n.r`, and both arms share `pmin` -- if that cap bound, the two numbers
would be identical by construction and the comparison would prove nothing. It
does NOT bind at any depth tested (reported p is 1.5e-4 against a cap of
9.4e-4; 1.1e-14 against 7.4e-14; 8.7e-31 against 2.0e-30). The agreement is
real.

Limits, stated: 6 unique multi-rho cells from one fixture; the deep points are
manufactured by scaling; and below ~1e-30 the `pchisq` term underflows to zero
so the integral carries 100% of the p, yet the arms still agree to 5 significant
figures -- **observed and NOT explained**. That regime is far outside anything
this pipeline reports, so it is recorded rather than chased.

## 2. The fastTest escalation lever: 1.054x, not the 19.4% it looks like

The profile run makes 249 `mainRegionInCPP` calls for 210 regions: **39
escalations, 10.86 s of the 52.44 s of C++ (19.4%)**, and escalated regions
cost 1.30x an average first pass (they are the signal-bearing ones). That
looked like the largest remaining lever, and it is bit-exact by construction --
reuse, not recomputation.

**It is not, once the bucket is split.** The escalation re-tests with the
NON-fast path, so its tests are computing something different and legitimately
cannot be reused:

| escalation work | seconds |
|---|---|
| tests (`test` + `urTest`), legitimately re-run fast -> non-fast | 5.87 |
| genotype prep (everything else), the REUSABLE part | 4.99 |

Reusing the prep only: **96.69 -> 91.70 s = 1.054x at N = 50,000**, and
**1.095x at N = 400,000** (the prep is O(N), so it scales). For a change to
the most intricate control flow in the package.

## Standing conclusion for the SKAT-O path

Every lever is now costed against a measured profile:

| lever | N=50k | N=400k |
|---|---|---|
| whole R kernel free (unreachable ceiling) | 1.72x | 1.09x |
| all SKAT-O-only C++ free (unreachable ceiling) | 1.30x | 1.45x |
| escalation prep reuse | 1.054x | 1.095x |
| carrier-space Gram | 1.034x | 1.062x |
| assertions/instrumentation, all of it | ~1.01x | -- |

**#158 took the available win (2.14x, measured). Nothing large is left.** Three
times now a headline share has deflated on contact with the actual bucket
(carrier-space Phi 54% -> 3.5%, escalation 19.4% -> 9.5%, armadillo bounds
checks assumed absent and in fact present but ~1%). Treat "there is a big win
hiding in the O(N) buckets" as refuted, not unexplored.
