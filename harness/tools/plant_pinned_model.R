#!/usr/bin/env Rscript
# Plant the fitted values of a separated group at a chosen probability
# (LEDGER #174). Usage:
#   Rscript plant_pinned_model.R <in.rda> <out.rda> <mu>
#
# The separated group of the pinned fixture is fitted to ~1e-9..1e-7 because
# that is where R's glm stops iterating (harness/tools/make_pinned_fixture.py),
# not because the data say so: one more iteration moves every one of them by a
# factor of e. A fit that kept going would put them wherever it was stopped,
# down to R's clamp at machine epsilon -- one All of Us model has fitted values
# there. This writes the model that fit would have produced: the group's fitted
# values at <mu>, and everything step 2 reads that is a function of them
# (residuals, linear predictors, the projection matrices of obj.noK) recomputed
# by the package's own constructor. theta, the coefficients and the variance
# ratios are left alone: the group carries no information at either value.
#
# The group is every sample whose fitted value is at or below 1e-6 and whose
# outcome is 0; the script refuses a model in which that is nobody.
suppressMessages(library(SAIGE))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
value <- as.numeric(args[3])
stopifnot(is.finite(value), value > 0, value < 1e-6)
e <- new.env()
load(args[1], envir = e)
modglmm <- e$modglmm
mu <- as.numeric(modglmm$fitted.values)
y <- as.numeric(modglmm$y)
sep <- mu <= 1e-6 & y == 0
if (!any(sep)) stop("no separated sample (fitted value <= 1e-6, outcome 0) in ", args[1])
mu[sep] <- value
modglmm$fitted.values[] <- mu
modglmm$residuals[] <- y - mu
modglmm$linear.predictors[sep] <- stats::qlogis(value)
modglmm$obj.noK <- SAIGE:::ScoreTest_NULL_Model(mu, mu * (1 - mu), y, modglmm$X)
save(modglmm, file = args[2])
cat("planted", sum(sep), "fitted values at", format(value), "in", args[2], "\n")
