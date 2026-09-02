# Standalone probes

Small, self-contained experiments whose results are quoted in
`docs/design/step2-threading-verdict.md`. They are not part of any gate; they
exist so a claim in the docs can be re-checked rather than believed.

## `lchoose_offthread.cpp`

Asks whether `Rf_lchoose` -- the dependency that makes the ER exact test
main-thread-only, and so the reason in-process threading defers every binary
marker with `MAC <= max_MAC_use_ER` -- is genuinely unsafe off the main thread,
or merely tripping R's C-stack check.

```sh
cd harness/probes
PKG_CXXFLAGS="-std=c++14" pixi run R CMD SHLIB lchoose_offthread.cpp
pixi run Rscript -e 'dyn.load("lchoose_offthread.so");
  print(matrix(.Call("probe_lchoose", FALSE), ncol=2, byrow=TRUE))'   # aborts
pixi run Rscript -e 'dyn.load("lchoose_offthread.so");
  print(matrix(.Call("probe_lchoose", TRUE),  ncol=2, byrow=TRUE))'   # works
```

**Result (2026-09-02).** With the check active the worker call aborts the
session -- `C stack usage 17556326859020 is too close to the limit`, reproducing
S6's finding, the huge number being `main_sp - worker_sp`. With
`R_CStackLimit = (uintptr_t)-1` the same call **succeeds and returns bits
identical to the main thread** for (1,0), (1,1), (5,2) and (999,98).

So the restriction is R's stack check, not `lchoose`'s arithmetic. What that
does and does not license is in the verdict doc; in short it makes the ER
deferral look removable, but it does not on its own prove the rest of the ER
machinery is thread-safe.

## `offthread_localise.cpp` — WHICH call is the blocker

One entry point per candidate, each run in its own process (an R error raised on
a non-main thread is not something to recover from and keep using the session).
`drive_localise.R` drives it. **Include the no-op entry**: the first version of
this experiment reported "aborted" for everything including the empty worker,
which was a broken shell quoting layer, not a finding. A negative control is the
only reason that was caught.

**Result (2026-09-02), stack check ON:**

| off-thread call | result |
|---|---|
| empty worker (control) | OK |
| `Rf_lgammafn` | OK |
| `Rf_pnorm5` | OK |
| `Rf_lbeta` | **OK** |
| `Rf_lchoose(1,0)` | ABORT |
| `Rf_lchoose(5,2)` | ABORT |
| `Rf_choose(5,2)` | ABORT |

The decisive row is `lchoose(1,0)`: it returns 0 through an early return that
makes NO math-library call, and still aborts. So the trigger is at the ENTRY of
`choose`/`lchoose` -- the two functions in this list that recurse -- and not in
the arithmetic. `lbeta`, which is what `lchoose` actually computes with, is fine.

## `lchoose_replacement.cpp` — the fix, and its proof

`HyperGeo::lCombinations` calls `lchoose(int, int)` with `k > n` already
short-circuited, so the reachable domain is integer `n >= 0`, integer `k <= n`.
On that domain R's `choose.c` reduces to early returns plus
`lfastchoose(n,k) = -log(n+1) - lbeta(n-k+1, k+1)`, with its two recursive arms
unrollable (`n < 0` cannot occur; the symmetry arm lands on `n-k` in {0,1},
which are the `k < 2` early returns).

**Result (2026-09-02).** `drive_exhaustive.R`: over n in [0,2000], k in [-3,n] --
**2,009,004 pairs, ZERO bit-differences** against `Rf_lchoose` (compared with
`memcmp`, so NaN/Inf/-0.0 are held to the same standard). `drive_offthread.R`:
the same 2,009,004 values computed on a WORKER thread **with R's stack check
left ON** -- no abort, zero bit-differences against the main thread.

So the ER deferral can be removed **without** touching `R_CStackLimit`, i.e.
without giving up R's genuine stack-overflow protection, and byte-identical
goldens follow by construction.

**What is still open**: this removes the `lchoose` dependency only. `SL_setseed`
remains on the resampling arms (S2 made it conditional, and MAC <= 4 markers do
not take them), and the rest of the ER machinery -- `SL_calloc`, the `Binary_*`
statics -- has never been TSAN-tested with the deferral disabled. That test is
the safety case; `harness/tsan.sh` is the instrument.
