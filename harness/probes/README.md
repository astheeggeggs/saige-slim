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
