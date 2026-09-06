# Several null models per region (step 4, 2026-09-06), end to end through the
# CLI. The claim is the marker path's, now for regions: EACH model's region
# file, single-in-group file and marker list are byte-identical to what that
# model's own single-model run produces. One code path -- prepare each marker
# once, test it for every model, one region call per gene -- checked against
# itself with one model; not two code paths compared (the standing rule
# against which still holds). A one-row manifest writes the unsuffixed file
# and must reproduce the P2 / P5 goldens byte for byte (negative control for
# the whole feature). Models on different sample sets are refused before any
# state reaches C++.
#
# Tier: standard (runs the CLI).

RM_G <- c("--bedFile=./input/genotype_100markers.bed",
          "--bimFile=./input/genotype_100markers.bim",
          "--famFile=./input/genotype_100markers.fam",
          "--AlleleOrder=alt-first")
RM_SG <- c(
  paste0("--sparseGRMFile=./output/sparseGRM_relatednessCutoff_0.125_1000_",
         "randomMarkersUsed.sparseGRM.mtx"),
  paste0("--sparseGRMSampleIDFile=./output/sparseGRM_relatednessCutoff_0.125_",
         "1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt")
)
RM_GROUP <- c("--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
              "--groupFile=./input/group_new_snpid.txt",
              "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
              "--maxMAF_in_groupTest=0.0001,0.001,0.01",
              "--is_output_markerList_in_groupTest=TRUE")

# top-level aliases: lintr cannot see helper-saige.R's objects from inside a
# function body (the same device as test-step2-region-threads.R)
rm_cli <- run_saige_cli
RM_STEP1 <- file.path(saige_goldens, "step1")
RM_STEP2 <- file.path(saige_goldens, "step2")

rm_model <- function(prefix) {
  c(rda = file.path(RM_STEP1, paste0(prefix, ".rda")),
    vr = file.path(RM_STEP1, paste0(prefix, ".varianceRatio.txt")))
}

# name -> step-1 golden prefix (or a named c(rda=, vr=) pair)
rm_manifest <- function(entries) {
  path <- tempfile("region_manifest_", fileext = ".txt")
  lines <- vapply(names(entries), function(nm) {
    p <- entries[[nm]]
    if (is.character(p) && length(p) == 1) p <- rm_model(p)
    paste(nm, p[["rda"]], p[["vr"]])
  }, character(1))
  writeLines(unname(lines), path)
  path
}

rm_out <- function(tag) {
  out <- file.path(tempfile(paste0("step2_rm_", tag, "_")), "out.txt")
  dir.create(dirname(out), recursive = TRUE)
  out
}

rm_run <- function(args, out) {
  log <- tempfile(fileext = ".log")
  st <- rm_cli("step2_SPAtests.R", c(args, paste0("--SAIGEOutputFile=", out)), log = log)
  # keep the "log" attribute as the PATH expect_cli_ok prints on failure; the
  # lines go under their own name
  structure(st, loglines = readLines(log, warn = FALSE))
}

# every file a region run writes for one output prefix, present or not
rm_files <- function(prefix) {
  paste0(prefix, c("", ".singleAssoc.txt", ".markerList.txt"))
}

expect_same_region_files <- function(a, b) {
  fa <- rm_files(a)
  fb <- rm_files(b)
  for (j in seq_along(fa)) {
    testthat::expect_identical(file.exists(fa[j]), file.exists(fb[j]), info = fa[j])
    if (file.exists(fa[j])) {
      testthat::expect_identical(readBin(fa[j], "raw", file.info(fa[j])$size),
                                 readBin(fb[j], "raw", file.info(fb[j])$size),
                                 info = fa[j])
    }
  }
}

# The two step-1 goldens fitted on the SAME 1,000 samples (the marker path's
# multi-model test uses the same pair); the quantitative model's VR file has
# only `null` rows, so these runs pass no sparse GRM (as P36/P37 do).
RM_TWO <- list(binary = "binary_sparseGRM_a", quant = "quant_sparseGRM")

rm_single_args <- function(prefix, extra) {
  p <- rm_model(prefix)
  c(RM_G, RM_GROUP, paste0("--GMMATmodelFile=", p[["rda"]]),
    paste0("--varianceRatioFile=", p[["vr"]]), extra)
}

test_that("two models per region (SKAT-O, fastTest escalation) reproduce their single-model runs", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  extra <- c("--is_fastTest=TRUE", "--is_single_in_groupTest=TRUE")
  man <- rm_manifest(RM_TWO)
  out <- rm_out("skato")
  st <- rm_run(c(RM_G, RM_GROUP, extra,
                 paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man)), out)
  expect_cli_ok(st)
  testthat::expect_true(any(grepl("Region tests will run 2 null model(s) over one sample set",
                                  attr(st, "loglines"), fixed = TRUE)))
  for (nm in names(RM_TWO)) {
    ref <- rm_out(paste0("skato_", nm))
    expect_cli_ok(rm_run(rm_single_args(RM_TWO[[nm]], extra), ref))
    expect_same_region_files(paste0(out, "_", nm), ref)
    testthat::expect_true(file.exists(paste0(out, "_", nm, ".index")))
  }
})

# Quantitative model FIRST: an escalation call whose subset starts at the
# binary model then has a different first trait type from the run's, which is
# what the five-trait sim50k manifest tripped on (2026-09-06) -- the trait type
# handed to C++ must be the subset's own.
test_that("two models per region, quantitative first (SKAT-O, fastTest) reproduce single runs", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  extra <- c("--is_fastTest=TRUE", "--is_single_in_groupTest=TRUE")
  two <- RM_TWO[c("quant", "binary")]
  man <- rm_manifest(two)
  out <- rm_out("skato_qfirst")
  st <- rm_run(c(RM_G, RM_GROUP, extra,
                 paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man)), out)
  expect_cli_ok(st)
  for (nm in names(two)) {
    ref <- rm_out(paste0("skato_qfirst_", nm))
    expect_cli_ok(rm_run(rm_single_args(two[[nm]], extra), ref))
    expect_same_region_files(paste0(out, "_", nm), ref)
  }
})

test_that("two models per region (BURDEN, fastTest, 4 threads) reproduce their single-model runs", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  extra <- c("--r.corr=1", "--is_fastTest=TRUE", "--nThreads=4")
  man <- rm_manifest(RM_TWO)
  out <- rm_out("burden")
  st <- rm_run(c(RM_G, RM_GROUP, extra,
                 paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man)), out)
  expect_cli_ok(st)
  testthat::expect_true(any(grepl("Region marker loop will use 4 in-process threads",
                                  attr(st, "loglines"), fixed = TRUE)))
  for (nm in names(RM_TWO)) {
    ref <- rm_out(paste0("burden_", nm))
    expect_cli_ok(rm_run(rm_single_args(RM_TWO[[nm]], extra), ref))
    expect_same_region_files(paste0(out, "_", nm), ref)
  }
})

test_that("a one-row manifest writes the unsuffixed file and reproduces the P2 and P5 goldens", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  one <- rm_manifest(list(only = c(rda = "./output/example_binary_sparseGRM.rda",
                                   vr = "./output/example_binary_sparseGRM.varianceRatio.txt")))
  # P2: SKAT-O with fastTest escalation
  out <- rm_out("p2")
  st <- rm_run(c(RM_G, RM_SG, RM_GROUP[!grepl("markerList", RM_GROUP)], "--is_fastTest=TRUE",
                 paste0("--GMMATmodel_varianceRatio_multiTraits_File=", one)), out)
  expect_cli_ok(st)
  testthat::expect_false(file.exists(paste0(out, "_only")))
  expect_output_matches_golden(file.path(RM_STEP2, "plink_group_fastTest.txt"), out)
  # P5: BURDEN
  out <- rm_out("p5")
  st <- rm_run(c(RM_G, RM_SG, RM_GROUP[!grepl("markerList", RM_GROUP)],
                 "--r.corr=1", "--is_fastTest=FALSE",
                 paste0("--GMMATmodel_varianceRatio_multiTraits_File=", one)), out)
  expect_cli_ok(st)
  expect_output_matches_golden(file.path(RM_STEP2, "plink_group_burden.txt"), out)
})

test_that("models on different sample sets are refused on a region run, before C++", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # No two committed models disagree about sampleID (every step-1 golden and
  # the bundled example model carry the same 1,000-sample vector), so one is
  # manufactured exactly as test-step2-multimodel.R does: binary_sparseGRM_a
  # with its last 10 samples dropped, obj.noK rebuilt for the subset.
  sub <- local({
    e <- new.env()
    load(file.path(RM_STEP1, "binary_sparseGRM_a.rda"), envir = e)
    m <- e$modglmm
    n <- length(m$sampleID)
    keep <- seq_len(n - 10L)
    for (nm in names(m)) {
      x <- m[[nm]]
      if (is.matrix(x)) {
        if (nrow(x) == n) m[[nm]] <- x[keep, , drop = FALSE]
      } else if (is.atomic(x) && length(x) == n) {
        m[[nm]] <- x[keep]
      }
    }
    mu <- as.vector(m$fitted.values)
    mu2 <- mu * (1 - mu)
    m$obj.noK <- SAIGE:::ScoreTest_NULL_Model(mu, mu2, m$y, m$X)
    modglmm <- m
    f <- tempfile(fileext = ".rda")
    save(modglmm, file = f)
    f
  })
  man <- rm_manifest(list(binsub = c(rda = sub, vr = rm_model("binary_sparseGRM_a")[["vr"]]),
                          quant = "quant_sparseGRM"))
  out <- rm_out("refuse")
  st <- rm_run(c(RM_G, RM_GROUP, "--r.corr=1",
                 paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man)), out)
  testthat::expect_false(identical(as.integer(st), 0L))
  testthat::expect_true(any(grepl("SAME sample set", attr(st, "loglines"), fixed = TRUE)))
  testthat::expect_false(file.exists(paste0(out, "_binsub")))
})

# ---- deterministic escalation subsets -------------------------------------
# The two R defects sim50k caught on 2026-09-06 both needed an escalation
# subset that is a strict, non-leading part of the manifest: a list element
# deleted by `x[[k]] <- NULL` shifted the models after it, and the trait type
# handed to C++ was the run's first model's rather than the subset's. The
# bundled fixture cannot produce such a subset by chance, but
# --pval_cutoff_for_fastTest decides escalation per model from the fast pass's
# Cauchy p, which the single-model runs print: a cutoff BETWEEN the two
# models' p-values for GENE1 escalates exactly the model with the smaller p,
# and the manifest order then puts that model first or second. The log line
# "Non-fast test is performed" is printed once per region that escalates, so
# the subset is asserted, not assumed. Seconds, not the hour sim50k takes.
rm_cauchy_p <- function(out, region) {
  d <- utils::read.table(out, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  as.numeric(d$Pvalue[d$Region == region & d$Group == "Cauchy"])
}

rm_n_escalated <- function(st) {
  sum(grepl("Non-fast test is performed", attr(st, "loglines"), fixed = TRUE))
}

test_that("an escalation subset that is not the leading model reproduces the single-model runs", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # fast pass only (a cutoff of 0 would switch fastTest OFF, so the smallest
  # positive double): each model's fast-pass Cauchy p for GENE1
  probe <- c("--is_fastTest=TRUE", "--is_single_in_groupTest=TRUE",
             "--pval_cutoff_for_fastTest=1e-300")
  p <- vapply(names(RM_TWO), function(nm) {
    out <- rm_out(paste0("probe_", nm))
    st <- rm_run(rm_single_args(RM_TWO[[nm]], probe), out)
    expect_cli_ok(st)
    testthat::expect_identical(rm_n_escalated(st), 0L)
    rm_cauchy_p(out, "GENE1")
  }, numeric(1))
  testthat::expect_true(all(is.finite(p)))
  testthat::expect_true(p[["binary"]] != p[["quant"]])
  cutoff <- (p[["binary"]] + p[["quant"]]) / 2
  smaller <- names(p)[which.min(p)]
  larger <- setdiff(names(p), smaller)
  extra <- c("--is_fastTest=TRUE", "--is_single_in_groupTest=TRUE",
             paste0("--pval_cutoff_for_fastTest=", format(cutoff, digits = 17)))
  # the single runs at that cutoff: the smaller-p model escalates GENE1, the
  # other escalates nothing
  refs <- list()
  for (nm in names(RM_TWO)) {
    refs[[nm]] <- rm_out(paste0("esc_ref_", nm))
    st <- rm_run(rm_single_args(RM_TWO[[nm]], extra), refs[[nm]])
    expect_cli_ok(st)
    testthat::expect_identical(rm_n_escalated(st), if (nm == smaller) 1L else 0L, info = nm)
  }
  # the escalating model SECOND (the subset starts past the run's first model)
  # and FIRST (the subset is model 1 while model 2 keeps its fast pass)
  for (order in list(c(larger, smaller), c(smaller, larger))) {
    man <- rm_manifest(RM_TWO[order])
    out <- rm_out(paste0("esc_", paste(order, collapse = "_")))
    st <- rm_run(c(RM_G, RM_GROUP, extra,
                   paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man)), out)
    expect_cli_ok(st)
    testthat::expect_identical(rm_n_escalated(st), 1L, info = paste(order, collapse = ","))
    for (nm in order) {
      expect_same_region_files(paste0(out, "_", nm), refs[[nm]])
    }
  }
})
