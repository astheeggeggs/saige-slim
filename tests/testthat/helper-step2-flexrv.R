# Shared scaffolding for the FlexRV step-2 tests (test-step2-flexrv*.R).
# Nothing here runs at source time: fx_run() calls skip_unless_level("standard")
# before it spawns anything, so the fast tier pays only the definitions.
#
# The fixture, extdata/input/group_new_snpid_flexrv.txt, is
# group_new_snpid.txt (GENE1: 30 lof + 20 missense, GENE2: 20 lof + 30
# missense) plus one `score:REVEL` line per gene: LoF markers scored 1 (the
# paper's rule), missense markers carrying the REVEL values of
# group_new_snpid_multiWeights.txt, interior values in (0, 1) so every
# nonlinear score transform is live.
#
# The fixture ARITHMETIC (2026-09-09, measured with plink2 --freq on the model's
# 1,000 samples -- NOT the fam file's 10,000; step 2 subsets to the model):
# 2N = 2,000, so max MAF 0.01 admits MAC <= 20. At that cell GENE1 holds 21
# polymorphic markers with MAC 1..20 and GENE2 18 with MAC 1..13; at the
# default collapse cutoff (10) the real markers are rs4/rs17/rs31/rs33 (GENE1)
# and rs60/rs70/rs91 (GENE2), everything else is ultra-rare. Consequences:
#   - MAC25..MAC500 sit at or above the ingest bound 2*N*maxMAF = 20 and are
#     DROPPED before composition (LEDGER #142): a default run tests 112 sets, not 192;
#   - MAC1/3/5/10 select markers that are ALL collapsed at cutoff 10, so those
#     sets test the pseudo-marker alone (design doc, "What enters the test");
#   - with --MACCutoff_to_CollapseUltraRare=0 the whole MAC ladder is live on
#     real markers, which is the paper's unbinned quantitative run.
# 0.001 (MAC <= 2) would leave NO real marker at cutoff 10, which is why the
# tests and goldens pin --flexRV_maxMAF=0.01.

FX_G <- c("--bedFile=./input/genotype_100markers.bed",
          "--bimFile=./input/genotype_100markers.bim",
          "--famFile=./input/genotype_100markers.fam",
          "--AlleleOrder=alt-first")
FX_SG <- c(
  paste0("--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
         "randomMarkersUsed.sparseGRM.mtx"),
  paste0("--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_",
         "1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
)
FX_BIN <- c("--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
            "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt")
# the quantitative model runs WITHOUT the sparse GRM: its variance-ratio file
# has only `null` rows. A function, not a value: helper files are sourced in
# name order and helper-saige.R (which defines saige_goldens) comes after this
# one.
fx_quant <- function() {
  gold <- get("saige_goldens", envir = parent.frame())
  c(paste0("--GMMATmodelFile=", file.path(gold, "step1", "quant_sparseGRM.rda")),
    paste0("--varianceRatioFile=", file.path(gold, "step1", "quant_sparseGRM.varianceRatio.txt")))
}
FX_COMMON <- c("--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
               "--groupFile=./input/group_new_snpid_flexrv.txt")
# the FlexRV run proper: one annotation group, one max MAF, burden
FX_FLEX <- c("--flexRV_score=REVEL", "--flexRV_maxMAF=0.01",
             "--annotation_in_groupTest=lof:missense")
FX_BURDEN <- c("--r.corr=1", "--is_fastTest=FALSE")

# helper-saige.R's functions, bound here at top level: lintr cannot see them
# from inside a function body (the same aliasing test-step2-region-threads.R
# uses). This file is named so that it sorts AFTER helper-saige.R: testthat
# sources helpers in name order, and a helper-flexrv.R failed to load.
fx_skip <- skip_unless_level
fx_cli <- run_saige_cli
fx_read <- read_saige_table
fx_cli_ok <- expect_cli_ok

# One CLI run per distinct argument vector, cached per worker process.
fx_cache <- new.env(parent = emptyenv())
fx_run <- function(args, tag = "fx") {
  fx_skip("standard")
  key <- paste(args, collapse = " ")
  if (!is.null(fx_cache[[key]])) return(fx_cache[[key]])
  out <- file.path(tempfile(paste0("step2_", tag, "_")), "region.txt")
  dir.create(dirname(out), recursive = TRUE)
  log <- tempfile(fileext = ".log")
  # system2 goes through a shell: a transform spec's ';' would end the command
  args <- ifelse(grepl(";", args, fixed = TRUE), shQuote(args), args)
  status <- fx_cli("step2_SPAtests.R", c(args, paste0("--SAIGEOutputFile=", out)), log = log)
  region <- NULL
  if (identical(as.integer(status), 0L) && file.exists(out)) {
    region <- fx_read(out)
  }
  res <- list(status = status, log = readLines(log, warn = FALSE), path = out, region = region)
  fx_cache[[key]] <- res
  res
}
fx_ok <- function(res) {
  fx_cli_ok(res$status)
}

# The default binary BURDEN run (golden P58's command; 112 sets at N = 1,000).
fx_burden_bin <- function() {
  fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN), "bin")
}
# The paper's quantitative run: no ultra-rare collapse (golden P59's command).
fx_burden_quant <- function() {
  fx_run(c(FX_G, fx_quant(), FX_COMMON, FX_FLEX, FX_BURDEN, "--MACCutoff_to_CollapseUltraRare=0"),
         "quant")
}

# The rows of a region table, the pooled Cauchy row, and a log predicate.
fx_rows <- function(res) res$region[res$region$Group != "Cauchy", ]
fx_cauchy <- function(res) res$region[res$region$Group == "Cauchy", ]
fx_log_has <- function(res, pattern) any(grepl(pattern, res$log, fixed = TRUE))

# Numeric columns of two region tables compared cell for cell at a relative
# tolerance, rows matched by (Region, Group, max_MAF, key) where `key` names
# the row within the cell (a Weight label, or NULL for the legacy schema).
# The conditional BURDEN columns are included when both tables carry them.
#
# On the tolerance: the tables are PRINTED at 6-7 significant digits, so two
# runs that differ by one ULP can print one digit apart in the last place, and
# a cross-path comparison cannot honestly demand more than ~1e-5 relative
# (the numeric compare mode's bound). 1e-10 is only for two runs of ONE path.
fx_num_cols <- c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden", "SE_Burden",
                 "Pvalue_Burden_c", "Beta_Burden_c", "seBeta_Burden_c")
fx_expect_rows_equal <- function(a, b, cols = NULL, tol = 1e-10, label = "rows") {
  if (is.null(cols)) {
    cols <- intersect(fx_num_cols, intersect(colnames(a), colnames(b)))
  }
  testthat::expect_identical(nrow(a), nrow(b), label = paste(label, "row count"))
  for (cc in cols) {
    x <- suppressWarnings(as.numeric(a[[cc]]))
    y <- suppressWarnings(as.numeric(b[[cc]]))
    testthat::expect_identical(is.na(x), is.na(y), label = paste(label, cc, "NA pattern"))
    testthat::expect_equal(x[!is.na(x)], y[!is.na(y)], tolerance = tol, label = paste(label, cc))
  }
}
# The allele accounting is weight-independent (real markers) or set by the
# collapse (the pseudo-marker) and must agree to the printed digit.
fx_acct_cols <- c("MAC", "MAC_case", "MAC_control", "Number_rare", "Number_ultra_rare")
fx_expect_accounting_identical <- function(a, b, label = "accounting") {
  for (cc in intersect(fx_acct_cols, intersect(colnames(a), colnames(b)))) {
    testthat::expect_identical(a[[cc]], b[[cc]], label = paste(label, cc))
  }
}

# ---- the in-sample MAC per marker, from two independent sources ------------
# Every step-1 golden and the bundled example model carry the same 1,000-sample
# vector (test-step2-region-multimodel.R relies on it), and no bundled marker
# has a missing call, so MAC == AC_Allele2 and MAF == MAC / 2000 for every
# marker of every run here.
#  - fx_mac_from_golden(): P1's golden, --minMAC=1, every polymorphic marker.
#  - fx_mac_from_single(res): the run's own .singleAssoc.txt -- only markers
#    that survived the ultra-rare collapse appear there, so it is complete only
#    when the collapse cutoff is 0.
fx_mac_from_golden <- function() {
  gold <- get("saige_goldens", envir = parent.frame())
  d <- fx_read(file.path(gold, "step2", "plink_marker_binary.txt"))
  stopifnot(all(d$MissingRate == "0"))
  stats::setNames(as.numeric(d$AC_Allele2), d$MarkerID)
}
fx_mac_from_single <- function(res) {
  d <- fx_read(paste0(res$path, ".singleAssoc.txt"))
  d <- d[grepl("^rs", d$MarkerID), ]
  stopifnot(all(d$MissingRate == "0"))
  stats::setNames(as.numeric(d$AC_Allele2), d$MarkerID)
}

# ---- the external route: weights written out as custom `weight:NAME` lines --
# One group file with the fixture's var/anno lines and K weight lines per
# region, the weights computed in R by `weightFun(score, isLof, MAF, MAC)`
# (a q x K matrix) and written at 17 digits. A marker absent from `mac` is
# monomorphic in the model's samples (MAC 0): it fails QC in every run, so its
# weight never enters a test. Run as ordinary custom sets, this is the
# reference the internal transform expansion is reconciled against.
fx_weight_groupfile <- function(mac, setNames, weightFun,
                                path = tempfile("flexrv_external_", fileext = ".txt")) {
  ext <- get("saige_extdata", envir = parent.frame())
  fixture <- readLines(file.path(ext, "input", "group_new_snpid_flexrv.txt"))
  out <- character(0)
  for (i in seq(1, length(fixture), by = 3)) {
    var <- strsplit(fixture[i], " ")[[1]]
    anno <- strsplit(fixture[i + 1], " ")[[1]]
    sc <- strsplit(fixture[i + 2], " ")[[1]]
    stopifnot(var[2] == "var", anno[2] == "anno", sc[2] == "score:REVEL")
    ids <- var[-(1:2)]
    m <- unname(mac[ids])
    m[is.na(m)] <- 0
    W <- weightFun(as.numeric(sc[-(1:2)]), as.numeric(anno[-(1:2)] == "lof"), m / 2000, m)
    stopifnot(nrow(W) == length(ids), ncol(W) == length(setNames))
    out <- c(out, fixture[i], fixture[i + 1])
    for (k in seq_len(ncol(W))) {
      out <- c(out, paste(c(var[1], paste0("weight:", setNames[k]), format(W[, k], digits = 17)),
                          collapse = " "))
    }
  }
  writeLines(out, path)
  path
}
# The transform table's sets, under the internal run's own set names, so the
# two runs' Weight columns match label for label. flexrvWeightsR is the R
# route written from the paper's dbeta forms, not from the C++ closed forms.
fx_external_groupfile <- function(mac, scoreSpec = SAIGE::flexrvDefaultScoreSpec,
                                  mafSpec = SAIGE::flexrvDefaultMAFSpec, maxMAF = 0.01) {
  # deliberately WITHOUT N: the external route carries every set of the spec
  # as a custom weight line (192 by default), including the MAC rungs the
  # internal run drops at N = 1,000 (LEDGER #142) -- those lines are the second
  # computation that the dropped sets were exactly unw (fx_expect_int_eq_ext)
  tt <- SAIGE::flexrvTransformTable(scoreSpec, mafSpec, "REVEL", maxMAF = maxMAF)
  path <- fx_weight_groupfile(mac, tt$names, function(s, isLof, maf, m) {
    SAIGE::flexrvWeightsR(tt$table, s, isLof, maf, m, maxMAF)
  })
  list(path = path, names = tt$names)
}

# ---- run builders shared by the test-step2-flexrv-*.R files ----------------
# Defined HERE and not in the test files because lintr's object-usage check
# cannot see another file's objects from inside a function body (FX_G, fx_run
# and friends are visible in this file, where they are defined).

# The quantitative run with the collapse off (P59's command) plus
# --is_single_in_groupTest, so the run itself reports every marker's in-sample AC.
fx_quant_single <- function() {
  fx_run(c(FX_G, fx_quant(), FX_COMMON, FX_FLEX, FX_BURDEN, "--MACCutoff_to_CollapseUltraRare=0",
           "--is_single_in_groupTest=TRUE"), "quant-single")
}
# A plain (no --flexRV_score) custom-weight run of an external group file on
# the given model, with P58/P59's other flags.
fx_external_run <- function(path, model, extra, tag) {
  fx_run(c(FX_G, model, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
           paste0("--groupFile=", path), "--annotation_in_groupTest=lof:missense",
           "--maxMAF_in_groupTest=0.01", FX_BURDEN, extra), tag)
}
# The MAC ladder's second parameterization: a plain unweighted run whose
# --maxMAC_in_groupTest=t cells are the marker subsets the MAC<=t sets weight 1.
fx_maxmac_run <- function() {
  fx_run(c(FX_G, fx_quant(), FX_COMMON, "--annotation_in_groupTest=lof:missense",
           "--maxMAF_in_groupTest=0.01", "--maxMAC_in_groupTest=1,3,5,10",
           "--is_no_weight_in_groupTest=TRUE", "--MACCutoff_to_CollapseUltraRare=0", FX_BURDEN),
         "maxmac")
}
fx_ladder_row <- function(q, gene, t) {
  r <- fx_rows(q)
  r[r$Region == gene & r$Weight == paste0("REVEL_unw_MAC", t), ]
}
fx_cell_row <- function(plain, gene, t) {
  r <- fx_rows(plain)
  r[r$Region == gene & abs(as.numeric(r$max_MAF) - t / 2000) < 1e-12, ]
}
# the `null` variance ratios of the quantitative model, named by MAC category
fx_quant_vr <- function() {
  vr <- utils::read.table(sub("^--varianceRatioFile=", "", fx_quant()[2]), header = FALSE,
                          stringsAsFactors = FALSE)
  vr <- vr[vr$V2 == "null", ]
  stats::setNames(vr$V1, vr$V3)
}
# The reduced-spec SKAT-O run (golden P60's command): 2 x 2 sets at max MAF
# 0.02 with mafbeta rather than a MAC threshold, so no in-cell weight is 0
# (LEDGER #117).
fx_skato <- function(extra = character(0), tag = "skato-p60") {
  fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, "--flexRV_score=REVEL", "--flexRV_maxMAF=0.02",
           "--annotation_in_groupTest=lof:missense", "--flexRV_scoreTransforms=unw;sbeta:2",
           "--flexRV_MAFTransforms=unw;mafbeta:2", "--is_fastTest=FALSE", extra), tag)
}
fx_default_names <- function() {
  # the 112 sets a default run composes at N = 1,000 (LEDGER #142)
  SAIGE::flexrvTransformTable(SAIGE::flexrvDefaultScoreSpec, SAIGE::flexrvDefaultMAFSpec,
                              "REVEL", maxMAF = 0.01, N = 1000)$names
}
fx_fixture_lines <- function() {
  ext <- get("saige_extdata", envir = parent.frame())
  readLines(file.path(ext, "input", "group_new_snpid_flexrv.txt"))
}

# The whole-table reconciliation: an internal default-set run (112 at this N)
# against its external
# weight-line twin, rows at the printed precision, accounting exact, the
# Cauchy row, and a positive control that the comparison tells sets apart.
fx_expect_int_eq_ext <- function(int, ext, label) {
  a <- fx_rows(int)
  b <- fx_rows(ext)
  # the internal run composes 112 sets at N = 1,000 (LEDGER #142 drops
  # MAC25..MAC500); the external file carries all 192 as custom sets. The 112
  # shared names must agree row for row, and the 80 external rows the internal
  # run never composed must be their score transform's unw row -- the second
  # computation that what was dropped was a duplicate and nothing else.
  testthat::expect_identical(nrow(a), 224L, label = label)
  testthat::expect_identical(nrow(b), 384L, label = label)
  testthat::expect_identical(a$Weight, rep(fx_default_names(), 2), label = label)
  b2 <- b[b$Weight %in% a$Weight, , drop = FALSE]
  rownames(b2) <- NULL
  testthat::expect_identical(a$Weight, b2$Weight, label = label)
  fx_expect_rows_equal(a, b2, tol = 1e-5, label = label)
  fx_expect_accounting_identical(a, b2, label = label)
  # the internal pooled row is the CCT over the 112 shared rows (the external
  # run pools all 192 of its own, so its Cauchy row is a different statistic)
  for (g in c("GENE1", "GENE2")) {
    testthat::expect_equal(as.numeric(fx_cauchy(int)$Pvalue_Burden[fx_cauchy(int)$Region == g]),
                           get_CCT_pvalue(as.numeric(b2$Pvalue_Burden[b2$Region == g])),
                           tolerance = 1e-5, label = paste(label, "Cauchy", g))
  }
  strip <- function(d) {
    d$Weight <- NULL
    rownames(d) <- NULL
    d
  }
  dropped <- setdiff(unique(b$Weight), unique(a$Weight))
  testthat::expect_identical(length(dropped), 80L, label = label)
  # every dropped name is one of the five rungs (so the substitution below
  # never leaves a name unchanged and compares a row with itself)
  testthat::expect_true(all(grepl("_MAC(25|50|100|200|500)$", dropped)), label = label)
  for (w in dropped) {
    unw <- sub("_MAC(25|50|100|200|500)$", "_unw", w)
    testthat::expect_false(identical(unw, w), label = w)
    testthat::expect_identical(strip(b[b$Weight == w, ]), strip(b[b$Weight == unw, ]),
                               label = paste(label, w))
  }
  # unw_unw against the steepest score transform, and against the MAC10 mask
  p <- function(d, w) as.numeric(d$Pvalue_Burden[d$Weight == w])
  testthat::expect_false(isTRUE(all.equal(p(a, "REVEL_unw_unw"), p(b, "REVEL_sbeta20_unw"),
                                          tolerance = 1e-5)))
  testthat::expect_false(isTRUE(all.equal(p(a, "REVEL_unw_unw"), p(b, "REVEL_unw_MAC10"),
                                          tolerance = 1e-5)))
}

# ---- several null models in one pass (test-step2-flexrv-multimodel.R) ------
# The two step-1 goldens fitted on the SAME 1,000 samples; the quantitative
# model's VR file has only `null` rows, so no sparse GRM is passed.
FXM_MODELS <- c(binary = "binary_sparseGRM_a", quant = "quant_sparseGRM")

fxm_model <- function(prefix) {
  gold <- get("saige_goldens", envir = parent.frame())
  c(rda = file.path(gold, "step1", paste0(prefix, ".rda")),
    vr = file.path(gold, "step1", paste0(prefix, ".varianceRatio.txt")))
}
fxm_manifest <- function() {
  path <- tempfile("flexrv_manifest_", fileext = ".txt")
  writeLines(vapply(names(FXM_MODELS), function(nm) {
    p <- fxm_model(FXM_MODELS[[nm]])
    paste(nm, p[["rda"]], p[["vr"]])
  }, character(1)), path)
  path
}
fxm_single_model <- function(nm) {
  p <- fxm_model(FXM_MODELS[[nm]])
  c(paste0("--GMMATmodelFile=", p[["rda"]]), paste0("--varianceRatioFile=", p[["vr"]]))
}
fxm_expect_same_files <- function(a, b) {
  for (suffix in c("", ".singleAssoc.txt")) {
    fa <- paste0(a, suffix)
    fb <- paste0(b, suffix)
    testthat::expect_identical(file.exists(fa), file.exists(fb), label = fa)
    if (file.exists(fa)) {
      testthat::expect_identical(readLines(fa), readLines(fb), label = fa)
    }
  }
}
# `flags`: everything but the model. The manifest run's per-model files are
# <path>_<name>; each must be byte-identical to that model's own run's.
fxm_check <- function(flags, tag) {
  man <- fx_run(c(FX_G, FX_COMMON, flags,
                  paste0("--GMMATmodel_varianceRatio_multiTraits_File=", fxm_manifest())),
                paste0(tag, "-manifest"))
  fx_ok(man)
  testthat::expect_true(fx_log_has(man, "Region tests will run 2 null model(s) over one sample"))
  for (nm in names(FXM_MODELS)) {
    single <- fx_run(c(FX_G, FX_COMMON, flags, fxm_single_model(nm)), paste0(tag, "-", nm))
    fx_ok(single)
    fxm_expect_same_files(paste0(man$path, "_", nm), single$path)
    testthat::expect_true(file.exists(paste0(man$path, "_", nm, ".index")))
  }
  man
}

# ---- one transform set, K = 1 (test-step2-flexrv-k1.R) ---------------------
# quantitative, collapse off, so a comparison is between burdens over the
# same real markers with no pseudo-marker in play
fxk_flex <- function(extra, tag) {
  fx_run(c(FX_G, fx_quant(), FX_COMMON, FX_FLEX, "--MACCutoff_to_CollapseUltraRare=0", extra), tag)
}
fxk_plain <- function(extra, tag) {
  fx_run(c(FX_G, fx_quant(), FX_COMMON, "--annotation_in_groupTest=lof:missense",
           "--maxMAF_in_groupTest=0.01", "--MACCutoff_to_CollapseUltraRare=0", extra), tag)
}
