# LEDGER #158 regression guard on the FINISHED build (2026-09-19)

The item `next.md` carried as open -- "a null sweep per trait on the finished
build as the regression guard" -- is DONE. **It passes.** Build `0a70bb8b`,
clean `harness/build.sh`, 32 replicate nulls x {quantitative, binary K=0.01,
binary K=0.002} = 96 arms, 95 run here plus the pilot, every arm rc=0.

Data `harness/work/nullrep-finished-20260918/` (gitignored), summary
`guard_summary.rds`, scripts `dump_arm_finished.R` / `driver_finished.sh` /
`analyse_finished.R` beside it.

## What had to change in the harness, and why the old one would have lied

`harness/skato_pvalue/nullrep/dump_arm.R` patches
`SKAT:::Met_SKAT_Get_Pvalue`. **On the finished build that is no longer the
funnel every SKAT-O cell passes**: `get_SKAT_pvalue`
(`R/SAIGE_SPATest_Region_Func.R:411`) routes to `skatoScreenKernel()` whenever
the screen is on, and a screened cell never reaches SKAT's kernel. Re-running
that script as-is would have scored a SMALLER, BIASED cell set -- the
unscreened tail only -- and reported a pass.

So this guard patches `SAIGE:::get_SKAT_pvalue` instead, the one funnel
screened or not (`SAIGE_SPATest_Region.R`'s `skat_fun` default and
`skatoMemoKernel`'s `kernel` default both resolve it at call time, so the patch
is seen by every caller), and records per cell BOTH the shipped answer and
SKAT's own kernel on the SAME `(Score, Phi, r.corr)`. A miss is then counted
directly, not inferred. `nThreads = 1` keeps `canParallelSkato()` FALSE so no
fork carries the cells out of the process.

## Result

244,063 paired cells: 129,480 screened (`liu-screen`), 9,665 `davies`,
104,914 `single-rho`.

**Misses by this repo's definition (handoff-2026-09-17 section table: a miss is
`Davies p < threshold/10`, i.e. < 0.005 at the shipped threshold 0.05): ZERO.**
The 2026-09-18 claim reproduces on the finished build at 160x the cell count.

Finer grain, which the definition above hides and which nobody had measured at
this sample size: **149 of 129,480 screened cells (0.115%) have a true (Davies)
p below 0.05**, range [0.0475282, 0.0499698], reported as [0.0500048,
0.0527118], ratio shipped/upstream median 1.034 and max 1.063. That is entirely
inside Liu's known accuracy band and is the SAME phenomenon the 2026-09-17
table already recorded as "worst skipped Davies p 0.049 / 0.050" on ~1.5k
cells; at 244k cells the worst is 0.0475. Not a defect and not a miss -- the
consequence is reporting 0.0501 where the truth is 0.0498. **Zero of the 149
have a true p below 0.01**, so nothing is at stake anywhere a p-value is called.

Controls, both run before the number was believed:
- NEGATIVE: planting one miss on a correctly-routed cell gives 150 = real + 1.
  (The first version of this control planted on a cell that was ALREADY a miss
  and returned 149 -- it proved nothing. Plant on a non-miss.)
- PERMUTATION: pairing each shipped p against a permuted upstream p gives
  15,034, not 149, so the count is structured, not an artefact of the marginals.

## The null sweep cannot be the saddlepoint tail's guard, and this measures why

Of 244,063 cells: 58 below 1e-4, **3 unique below 1e-5**, ZERO below 1e-6. The
tail rows dump as 6 because of the fastTest escalation duplicate (all from
`y_nb002_r17`, cell idx 2274-2276 and 2288-2290, the same three cells twice).

Those 3 are the entire tail sample from 96 arms, and 2 of them take the
saddlepoint (q = 9 and q = 10; the third is `single-rho`, ratio 1.000). Their
shipped/Davies ratios are 1.254 and 1.118 -- **which says nothing about whether
the saddlepoint is right**, because Davies at SKAT's hardcoded `acc = 1e-6` is
not a valid reference below 1e-5; that is the premise #158 step 2 was built on.

Two things follow, and they point the same way:
1. **A null sweep structurally cannot validate the tail.** It cannot manufacture
   cells below 1e-5: 3 in 244,063.
2. **Both tail cells sit at q = 9 and q = 10**, i.e. inside the q = 4-10 band
   that carries 56% of this pipeline's Davies time and in which Kuonen's
   `o(n^-3/2)` error bound runs the WRONG way, with no published validation
   below n = 305 (`docs/state/research-2026-09-18-skato-speedups.md`, open
   question 2).
The exact-reference control that note proposes -- Ruben/Farebrother series or
extended-precision Imhof at small q, built in the LOWER tail to dodge the 1-F
cancellation -- is the only instrument that can settle it. This sweep is
evidence for doing it, not a substitute.

## Two facts about this run that are not results

- **The box slept mid-sweep.** `Low Power Sleep` at 22:57:43 on battery 1%,
  woken by AC attach at 01:40:09 (`pmset -g log`). Four arms were suspended for
  2h43m and resumed identically -- no correctness impact -- but their per-cell
  `t` / `t_upstream` fields (`proc.time()[["elapsed"]]`) have the sleep baked in.
  **Quote no timing from `y_nb002_r22`, `y_null_r23`, `y_nb01_r23`,
  `y_nb002_r23`.** `caffeinate -i -s -w <driver pid>` for the rest.
- This guard doubles the kernel work by construction (shipped answer plus
  SKAT's own on every cell), so its wall time is not a measurement of anything
  either.
