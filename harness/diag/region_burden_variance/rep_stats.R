# rep_stats.R <rep_dir> <prep.rds>: the burden score and its candidate null
# variances for the planted cell's burden variable, from the replicate's
# refitted null model -- pure R, no step 2: S = g~'(y - mu); V_naive =
# sum(W g~^2); V_sparse = g~' Sigma^-1 g~ (exact, symmetric); and the
# BURDEN-mode variance the C++ pseudo-variant test would use at this cell's
# MAC (> 20.5: null ratio, category 2, times V_naive).
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(TRUE)
prep <- readRDS(args[2])
load(file.path(args[1], "m.rda"))
stopifnot(identical(modglmm$sampleID, prep$sampleID))
u <- prep$burden
xv <- modglmm$obj.noK$XV
xxvx_inv <- modglmm$obj.noK$XXVX_inv
stopifnot(ncol(xv) == length(u), nrow(xxvx_inv) == length(u))
g_tilde <- u - as.numeric(xxvx_inv %*% (xv %*% u))
mu <- modglmm$fitted.values
w_mu <- mu * (1 - mu)
s <- sum(g_tilde * (modglmm$y - mu))
v_naive <- sum(w_mu * g_tilde^2)
sigma <- modglmm$theta[2] * prep$grm + Diagonal(x = 1 / w_mu)
v_sparse <- as.numeric(crossprod(g_tilde, solve(sigma, g_tilde)))
vr <- read.table(file.path(args[1], "m.varianceRatio.txt"), header = FALSE,
                 stringsAsFactors = FALSE)
vr_null2 <- vr$V1[vr$V2 == "null" & vr$V3 == 2]
vr_sparse2 <- vr$V1[vr$V2 == "sparse" & vr$V3 == 2]
saveRDS(list(S = s, V_naive = v_naive, V_sparse = v_sparse,
             V_burdenmode = vr_null2 * v_naive, V_sparse_vr = vr_sparse2 * v_sparse,
             tau2 = modglmm$theta[2], vr_null2 = vr_null2, vr_sparse2 = vr_sparse2,
             cases = sum(modglmm$y)),
        file.path(args[1], "rep_stats.rds"))
