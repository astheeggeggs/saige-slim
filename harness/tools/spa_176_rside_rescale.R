# LEDGER #176 review (2026-10-02): does the saddlepoint R runs itself,
# get_newPhi_scaleFactor_traitType, carry the same stretch as
# src/SAIGE_test.cpp:1050, and what does it do to the Phi rescale?
# Synthetic vectors only; needs the installed package. See
# docs/state/handoff-2026-10-02-176-review.md, evidence 5.
suppressMessages(library(SAIGE))
rescale <- SAIGE:::get_newPhi_scaleFactor_traitType
nSamples <- 5000

run <- function(secWeight, vr, label) {
  mu <- rep(0.01, nSamples) # a round SYNTHETIC value
  # raw burden vector: one heavy carrier, 81 light ones, non-carriers are 0
  g <- numeric(nSamples)
  g[1] <- 1
  g[2:82] <- secWeight / 81
  # the heavy carrier is a case, nobody else
  y <- numeric(nSamples)
  y[1] <- 1
  qSum <- sum(g * y)
  m1 <- sum(mu * g)
  var1 <- sum(mu * (1 - mu) * g^2)
  phi <- matrix(var1 * vr, 1, 1)
  score <- qSum - m1
  truth <- 0.01 # P(the atom fires), the exact upper tail here
  out <- rescale(qSum, mu, g, log(truth), score, phi, "SKAT-O", "binary")
  # the burden p the function used, by the same call it makes
  qAdj <- m1 + (qSum - m1) * sqrt(var1 / sum(phi))
  qInv <- 2 * m1 - qAdj
  pNoadj <- pchisq((qAdj - m1)^2 / var1, df = 1, lower.tail = FALSE, log.p = TRUE)
  conv <- TRUE
  pBurden <- SAIGE:::SPA_pval(mu, g, qAdj, qInv, pNoadj, .Machine$double.eps^0.25,
                              TRUE, "binary", conv)
  cat(sprintf(
    "%-34s sec weight %.3f VR %.3f | burden p used %.3e (exact tail %.3e) | Phi out/in %.4f\n",
    label, secWeight, vr, exp(pBurden), truth, out$val[1, 1] / phi[1, 1]
  ))
}

run(0.003, 1.00, "negative control (no stretch)")
run(0.003, 0.98, "one atom, stretched")
run(0.22, 0.98, "moderate secondary weight, stretched")
run(2.3, 1.00, "well spread, no stretch")
run(2.3, 0.98, "well spread, stretched")
