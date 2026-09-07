# planted_mle.R [n]: the maximum-likelihood burden coefficient per planted
# replicate, as the reference for the score estimator S/V. The phenotype and
# the true random effects are regenerated from the replicate's seed (sim_pheno.R,
# same RNG calls), so the conditional MLE (true b as an offset) is available,
# and the marginal one (b ignored) beside it.
args <- commandArgs(TRUE)
n_rep <- if (length(args) >= 1) as.integer(args[1]) else 100L
repo <- Sys.getenv("SAIGE_REPO", getwd())
out_dir <- Sys.getenv("DIAG_OUT", file.path(repo, "harness/work/diag-burden"))
prep <- readRDS(file.path(out_dir, "sim", "prep.rds"))
beta_b <- 0.05
n <- length(prep$eta_fixed)
res <- t(sapply(seq_len(n_rep), function(id) {
  set.seed(2000L + id)
  b <- as.numeric(prep$lmat %*% rnorm(n))
  eta <- prep$eta_fixed + b + beta_b * prep$burden
  y <- rbinom(n, 1, plogis(eta))
  cond <- glm(y ~ prep$burden, family = binomial, offset = prep$eta_fixed + b)
  marg <- glm(y ~ prep$burden, family = binomial, offset = prep$eta_fixed)
  c(cases = sum(y), mle_cond = unname(coef(cond)[2]), mle_marg = unname(coef(marg)[2]))
}))
# the regenerated phenotype must be the one the pipeline saw
log1 <- readLines(file.path(out_dir, "sim", "planted", "rep_0001", "log.txt"))
logged <- as.integer(sub(".*cases ([0-9]+).*", "\\1", grep("^seed", log1, value = TRUE)[1]))
stopifnot(identical(as.integer(res[1, "cases"]), logged))
cat(sprintf("planted arm, %d replicates, beta_b = %.3f\n", n_rep, beta_b))
for (cn in c("mle_cond", "mle_marg")) {
  cat(sprintf("  %-9s mean %.4f (%.3f x beta_b, SE %.3f)  sd %.4f\n", cn, mean(res[, cn]),
              mean(res[, cn]) / beta_b, sd(res[, cn]) / sqrt(n_rep) / beta_b, sd(res[, cn])))
}
