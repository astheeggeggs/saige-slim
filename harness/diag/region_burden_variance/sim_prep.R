# sim_prep.R <out.rds> [tau2]: everything the phenotype simulator needs, once.
# tau2 defaults to the model's own; a larger value amplifies the relatedness
# so the sparse-GRM and naive variances separate far enough to test.
# The null GLMM of the committed binary model: fixed part from its GLM fit,
# random effects b ~ N(0, tau2 * sparse GRM), and SAIGE's own weighted burden
# variable for the planted-effect arm (cell 4 of the R0 dump, GENE1
# missense;lof 0.01, residualised and beta-weighted as the region path forms it).
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(TRUE)
repo <- Sys.getenv("SAIGE_REPO", getwd())
out_dir <- Sys.getenv("DIAG_OUT", file.path(repo, "harness/work/diag-burden"))
load(file.path(repo, "extdata/output/example_binary_sparseGRM.rda"))
grm_name <- "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx"
grm_file <- file.path(repo, "extdata/output", grm_name)
grm <- readMM(grm_file)
ids <- readLines(paste0(grm_file, ".sampleIDs.txt"))
idx <- match(modglmm$sampleID, ids)
stopifnot(!anyNA(idx))
grm <- as(as(grm, "generalMatrix"), "CsparseMatrix")[idx, idx]
tau2 <- if (length(args) >= 2) as.numeric(args[2]) else modglmm$theta[2]
# covariates enter this model as an offset (isCovariateOffset), the only
# fitted fixed effect being the intercept
stopifnot(isTRUE(modglmm$isCovariateOffset))
eta_fixed <- as.numeric(modglmm$offset) + as.numeric(modglmm$coefficients)[1]
n <- length(modglmm$y)
stopifnot(length(eta_fixed) == n)
e <- eigen(tau2 * as.matrix(grm), symmetric = TRUE)
lmat <- e$vectors %*% diag(sqrt(pmax(e$values, 0)))
burden <- readRDS(file.path(out_dir, "R0", "phi_04.rds"))$g.sum
stopifnot(length(burden) == n)
pheno_name <- "pheno_1000samples.txt_withdosages_withBothTraitTypes.txt"
pheno_file <- file.path(repo, "extdata/input", pheno_name)
pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
stopifnot(setequal(pheno$IID, modglmm$sampleID))
cat(sprintf("n = %d  cases = %d  tau2 = %.4f  min eigen = %.3g  prevalence(fixed) = %.4f\n",
            n, sum(modglmm$y), tau2, min(e$values), mean(plogis(eta_fixed))))
saveRDS(list(sampleID = modglmm$sampleID, eta_fixed = eta_fixed, lmat = lmat,
             burden = burden, pheno = pheno, grm = grm, tau2 = tau2),
        args[1])
