# sim_pheno.R <prep.rds> <seed> <beta_b> <out_pheno.txt>
# y ~ Bernoulli(plogis(eta_fixed + b + beta_b * burden)), b ~ N(0, tau2 GRM).
args <- commandArgs(TRUE)
prep <- readRDS(args[1])
seed <- as.integer(args[2])
beta_b <- as.numeric(args[3])
set.seed(seed)
n <- length(prep$eta_fixed)
b <- as.numeric(prep$lmat %*% rnorm(n))
eta <- prep$eta_fixed + b + beta_b * prep$burden
y <- rbinom(n, 1, plogis(eta))
pheno <- prep$pheno
pheno$y_binary <- y[match(pheno$IID, prep$sampleID)]
write.table(pheno, args[4], quote = FALSE, row.names = FALSE, sep = " ")
cat(sprintf("seed %d  beta_b %.3f  cases %d  var(b) %.3f\n", seed, beta_b, sum(y), var(b)))
