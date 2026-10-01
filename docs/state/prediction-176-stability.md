# Prediction: is LEDGER #176 saddlepoint INSTABILITY on a concentrated vector? (2026-10-01)

Written BEFORE the experiment is run, so it can fail.

## First, the row's two symptoms are ONE symptom

`SE_Burden = |BETA_Burden / qnorm(p/2)|` -- it is BACK-DERIVED from the p-value
(`R/SAIGE_SPATest_Region_Func.R:490`), so it carries no information the p does
not. With `BETA_Burden` identical at 83.780 and 83.788:

    SE ratio  = 30.54 / 9.78  = 3.123
    |Z| ratio = 8.57  / 2.74  = 3.124      (the row's own 8.56 and 2.74)

**So "the 3.12 SE factor" and "|Z| is bimodal" are the same observation**, and
the row has exactly ONE anomaly, not two: two nearly identical weight sets got
very different p-values from the same statistic. Anything that explains the
p-values explains the SE.

## The hypothesis

The cell is CONCENTRATED (one atom carries the weight) and ABOVE the 13-carrier
cap, so #172's exact-test rescue cannot fire and the saddlepoint runs. If the
saddlepoint is UNSTABLE in that regime, two nearly identical vectors can land on
very different roots -- which is exactly what `sbeta10` against `sbeta20` is,
since `BETA_Burden` agreeing to four significant figures says the two vectors are
the same to ~1e-4 relative.

## The experiment, on the bundled fixture -- NO AoU data

`input/group_new_snpid_withWeights_abovecap.txt` (built 2026-10-01 for #176's
counter) already produces both regimes in one run:

- cells with `n_material` 15..94 and `n_eff` <= 4 -- CONCENTRATED, ABOVE the cap,
  saddlepoint. **The test.**
- cells with `n_material` 5 and 10 -- concentrated, BELOW the cap, routed to the
  exact test. **The built-in negative control**, in the same run, same genotypes.

Run a LADDER of group files whose only difference is the non-atom weight, through
+-10% of 0.01 with fine steps at +-0.1%, and read `Pvalue_Burden` per cell.

## PREDICTIONS

1. **If the hypothesis holds**: the above-cap cells' p-values move NON-smoothly
   across the ladder -- jumps of orders for a 0.1% weight change -- while the
   routed cells move smoothly. That reproduces #176 with no AoU data and
   localises it to the saddlepoint on concentrated vectors.
2. **If everything moves smoothly**: the hypothesis is FALSE. #176 is not
   saddlepoint instability, and the amr observation needs another account -- the
   first candidate then being that the two sets' vectors are NOT as similar as
   `BETA_Burden` agreeing to 4 s.f. suggests, because BETA is a ratio
   (`sum(Score)/sum(Phi)`) and is insensitive to a common scaling of the weights.
3. **Intermediate**: smooth but with a step at a particular weight -- that would
   point at a BRANCH being crossed (the `|Z| < 2` cutoff in
   `get_newPhi_scaleFactor_traitType`, which switches between the normal
   approximation and the saddlepoint), not at instability. Note the row's three
   |Z| values are 0.51, 8.56 and 2.74, and **0.51 is below that cutoff** -- so at
   least one of amr's three sets took the normal approximation and the other two
   took the saddlepoint. That alone could produce a bimodal |Z|, and it is the
   cheapest explanation on the table.

## What would NOT follow from outcome 1

That the amr gene's p is wrong in a particular direction. Instability says the
answer is unreliable, not which way. Establishing the TRUE p for a concentrated
above-cap cell needs a reference -- the exact test with the cap lifted, or a
Monte Carlo over the carriers -- which is a separate step.
