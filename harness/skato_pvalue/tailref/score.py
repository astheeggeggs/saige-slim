import sys, csv, numpy as np
sys.path.insert(0, "/private/tmp/claude-501/-Users-dpalmer-Repositories-SAIGE/683ae45d-f885-4566-ab5a-7c188958d73b/scratchpad/tailref")
from ref_ruben import ruben_survival

rows = list(csv.DictReader(open(sys.argv[1] + "/eval.csv")))
out = []
for r in rows:
    lam = np.array([float(x) for x in r["lambda"].split("|")])
    cond = float(r["cond"]); q = float(r["q"])
    if cond > 5000:
        S, R, ok, K = (np.nan, np.nan, False, 0)      # reference infeasible by design
    else:
        S, R, ok, K = ruben_survival(lam, q, Kmax=20000, rtol=1e-6)
    out.append(dict(r, p_ref=S, ref_bound=R, ref_ok=int(ok), ref_K=K))
w = csv.DictWriter(open(f"{sys.argv[1]}/scored.csv", "w"), fieldnames=list(out[0].keys()))
w.writeheader(); w.writerows(out)
n_ok = sum(o["ref_ok"] for o in out)
print(f"rows {len(out)}  reference converged (rtol 1e-6): {n_ok}  refused: {len(out)-n_ok}")
