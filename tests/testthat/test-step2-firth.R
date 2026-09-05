# The LIVE Firth path: --is_Firth_beta=TRUE routes qualifying markers through
# SAIGEClass::fast_logistf_fit_simple (src/SAIGE_test.cpp:761), the last arm of
# the getMarkerPval cascade (score -> SPA -> ER -> Firth). No golden run sets
# the flag, so before this file the entire arm executed zero times in the suite
# while getMarkerPval itself ran ~6000 times.
#
# Firth is DELIBERATELY gated behind pCutoffforFirth (default 0.01) because the
# penalised fit is expensive: it is meant to refine effect sizes for the few
# markers whose (potentially inflated) p-value is already small. These tests
# therefore run it at 0.05 -- a realistic cutoff that fires on 5 of the 70
# bundled markers. Do NOT raise the cutoff to force Firth onto every marker:
# it is the expensive path, and applying it to null markers exercises a regime
# the estimator was never intended for (see LEDGER #26).
#
# Note this is a DIFFERENT implementation from the Rcpp-exported free function
# fast_logistf_fit characterised in test-unit-firth.R; they share a name and a
# purpose, not code.
#
# Tier: standard (runs the CLI).

FIRTH_ARGS <- function(out, extra = character(0)) c(
  "--bedFile=./input/genotype_100markers.bed",
  "--bimFile=./input/genotype_100markers.bim",
  "--famFile=./input/genotype_100markers.fam",
  "--AlleleOrder=alt-first",
  "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
  paste0("--GMMATmodelFile=", file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
  paste0("--varianceRatioFile=",
         file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt")),
  "--is_output_moreDetails=TRUE",
  extra,
  paste0("--SAIGEOutputFile=", out)
)

firth_cache <- new.env(parent = emptyenv())
firth_run <- function(extra = character(0)) {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # paste(character(0), collapse = " ") is "", and env[[""]] is an error
  key <- if (length(extra)) paste(extra, collapse = " ") else "<none>"
  if (!is.null(firth_cache[[key]])) return(firth_cache[[key]])
  out <- file.path(tempfile("step2_firth_"), "marker.txt")
  dir.create(dirname(out), recursive = TRUE)
  expect_cli_ok(run_saige_cli("step2_SPAtests.R", FIRTH_ARGS(out, extra)))
  res <- read_saige_table(out)
  firth_cache[[key]] <- res
  res
}

pnum <- function(x) suppressWarnings(as.numeric(x))

# Independent, unadjusted log-odds-ratio from the genotype counts SAIGE itself
# reports -- a comparator that shares no code with the effect-size estimators.
unadj_logOR <- function(d) {
  ca <- 2 * pnum(d$N_case_hom) + pnum(d$N_case_het)
  ct <- 2 * pnum(d$N_case) - ca
  oa <- 2 * pnum(d$N_ctrl_hom) + pnum(d$N_ctrl_het)
  ot <- 2 * pnum(d$N_ctrl) - oa
  log(((ca + 0.5) / (ct + 0.5)) / ((oa + 0.5) / (ot + 0.5)))
}

marker_mac <- function(d) {
  N <- pnum(d$N_case) + pnum(d$N_ctrl)
  pmin(pnum(d$AC_Allele2), 2 * N - pnum(d$AC_Allele2))
}

test_that("is_Firth_beta at the default cutoff changes nothing on this fixture", {
  # The smallest p-value here is ~0.0100465, just above the 0.01 default, so no
  # marker qualifies and the flag is a no-op. This is the regression guard: if
  # the qualifying rule changes, the two runs stop matching.
  base <- firth_run()
  firth <- firth_run("--is_Firth_beta=TRUE")
  expect_identical(base, firth)
  expect_gt(min(pnum(base$p.value)), 0.01)
})

test_that("pCutoffforFirth selects exactly the markers at or below the cutoff", {
  base <- firth_run()
  firth <- firth_run(c("--is_Firth_beta=TRUE", "--pCutoffforFirth=0.05"))
  changed <- which(base$BETA != firth$BETA)
  qualifying <- which(pnum(base$p.value) <= 0.05)
  expect_identical(changed, qualifying)
  expect_gt(length(qualifying), 3)
  # everything above the cutoff is left exactly as the unpenalised run
  expect_identical(base$BETA[-qualifying], firth$BETA[-qualifying])
  expect_identical(base$SE[-qualifying], firth$SE[-qualifying])
})

test_that("Firth touches only the effect-size columns, never the score test", {
  base <- firth_run()
  firth <- firth_run(c("--is_Firth_beta=TRUE", "--pCutoffforFirth=0.05"))
  # the score test, allele accounting and sample counts are computed before the
  # Firth step and must be untouched by it
  for (col in c("MarkerID", "AC_Allele2", "AF_Allele2", "Tstat", "var",
                "p.value", "N_case", "N_ctrl", "N_case_het", "N_ctrl_het")) {
    expect_identical(base[[col]], firth[[col]], label = paste("unchanged:", col))
  }
  expect_false(identical(base$BETA, firth$BETA))
  expect_false(identical(base$SE, firth$SE))
})

test_that("in-regime Firth estimates agree with the score-based effect and the raw log-OR", {
  # For the well-powered qualifying markers the penalised fit must be
  # estimating the same quantity as the score-based BETA: same sign as both
  # that and the independent unadjusted log-OR, and the same order of
  # magnitude. Singletons are excluded -- at MAC <= 2 the score-based BETA is
  # itself unstable (one such marker reports 16.6 here), so neither estimator
  # is a usable reference.
  base <- firth_run()
  firth <- firth_run(c("--is_Firth_beta=TRUE", "--pCutoffforFirth=0.05"))
  lor <- unadj_logOR(base)
  mac <- marker_mac(base)
  i <- which(pnum(base$p.value) <= 0.05 & mac > 100)
  expect_gt(length(i), 2)

  b <- pnum(base$BETA)[i]
  f <- pnum(firth$BETA)[i]
  expect_identical(sign(f), sign(b))
  expect_identical(sign(f), sign(lor[i]))
  expect_true(all(abs(f / b) > 0.5 & abs(f / b) < 2))
  # SEs stay in the same regime too
  expect_true(all(pnum(firth$SE)[i] / pnum(base$SE)[i] > 0.5))
  expect_true(all(pnum(firth$SE)[i] / pnum(base$SE)[i] < 2))
})

# ---- LEDGER #26 FIXED: the cap is reported as NOT converged --------------------
# The one place the fixture reaches the drifting regime is by forcing the cutoff
# to 1, which the header above says never to do for a FIXTURE claim. It is done
# here on purpose, for the convergence FLAG only: fast_logistf_fit_simple used
# to set isfirthconverge = true when it hit maxit, so the run's
# "Firth approx was applied to N markers. M sucessfully converged." line counted
# the markers that had drifted at maxstep for 50 iterations (|BETA| = 750) as
# converged. Owner decision 2026-09-05: flag it. BETA and SE are unchanged.

firth_forced_log <- function() {
  skip_unless_level("standard")
  skip_if_no_goldens()
  out <- file.path(tempfile("step2_firth_cap_"), "marker.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- run_saige_cli("step2_SPAtests.R",
                      FIRTH_ARGS(out, c("--is_Firth_beta=TRUE", "--pCutoffforFirth=1")))
  expect_cli_ok(st)
  list(table = read_saige_table(out), log = readLines(attr(st, "log"), warn = FALSE))
}

test_that("LEDGER #26 FIXED: markers that hit the Firth iteration cap are not converged", {
  res <- firth_forced_log()
  line <- grep("^Firth approx was applied to ", res$log, value = TRUE)
  expect_length(line, 1L)
  pat <- "applied to ([0-9]+) markers\\. ([0-9]+) sucessfully converged"
  m <- regmatches(line, regexec(pat, line))[[1]]
  expect_length(m, 3L)
  applied <- as.integer(m[2])
  converged <- as.integer(m[3])
  expect_identical(applied, nrow(res$table))
  # the drifters: |BETA| at maxstep * maxit = 750 on this fixture's null markers
  drifted <- abs(pnum(res$table$BETA)) > 100
  expect_gt(sum(drifted), 0L)
  # every drifter hit the cap, so the converged count excludes at least them
  expect_lte(converged, applied - sum(drifted))
  expect_lt(converged, applied)
})
