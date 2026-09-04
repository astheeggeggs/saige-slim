# The multi-model marker scan, end to end through the CLI.
#
# The claim being gated is exact: EACH model's rows are byte-identical to what
# that model's own single-model run produces. That is not a claim about two code
# paths agreeing (the standing rule against which still holds); it is ONE code
# path -- prepare once, test each model -- checked against the goldens the same
# path produced when there was only ever one model to test.
#
# Tier: standard (runs the CLI). Its own file so its runs go in parallel with
# the other step-2 files rather than serialising behind them.

MM_PLINK <- c(
  "--bedFile=./input/genotype_100markers.bed",
  "--bimFile=./input/genotype_100markers.bim",
  "--famFile=./input/genotype_100markers.fam",
  "--AlleleOrder=alt-first"
)
MM_PGEN <- c(
  "--pgenFile=./input/genotype_100markers.pgen",
  "--pvarFile=./input/genotype_100markers.pvar",
  "--psamFile=./input/genotype_100markers.psam"
)
MM_COMMON <- c("--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
               "--is_output_moreDetails=TRUE")

mm_model_paths <- function(prefix) {
  c(
    file.path(saige_goldens, "step1", paste0(prefix, ".rda")),
    file.path(saige_goldens, "step1", paste0(prefix, ".varianceRatio.txt"))
  )
}

# name -> step-1 golden prefix
mm_manifest <- function(entries) {
  path <- tempfile(fileext = ".txt")
  lines <- vapply(names(entries), function(nm) {
    p <- mm_model_paths(entries[[nm]])
    paste(nm, p[1], p[2])
  }, character(1))
  writeLines(unname(lines), path)
  path
}

mm_run <- function(args, out) {
  st <- run_saige_cli("step2_SPAtests.R", c(args, paste0("--SAIGEOutputFile=", out)))
  expect_cli_ok(st)
  st
}

mm_log <- function(status) readLines(attr(status, "log"), warn = FALSE)

test_that("two models over one sample set reproduce their single-model goldens", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  man <- mm_manifest(list(binary = "binary_sparseGRM_a", quant = "quant_sparseGRM"))
  out <- file.path(tempfile("step2_mm_"), "multi.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- mm_run(c(MM_PLINK, MM_COMMON,
                 paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man)), out)

  # P1 and P32: the two goldens these exact flags produce one model at a time.
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "plink_marker_binary.txt"),
    paste0(out, "_binary")
  )
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "plink_marker_quantitative.txt"),
    paste0(out, "_quant")
  )

  # Both lines asserted for the LEDGER #30 reason: a manifest that failed to
  # arrive, or a scan that silently ran one model, would otherwise pass here.
  lg <- mm_log(st)
  expect_true(any(grepl("2 null model\\(s\\) in 1 sample-set group\\(s\\)", lg)))
  expect_true(any(grepl("Marker scan will test 2 models in 1 sample-set group", lg)))

  # Each model gets its OWN progress index, carrying its own byte count
  # (LEDGER #47): one shared count would be wrong for both.
  for (sfx in c("_binary", "_quant")) {
    idx <- readLines(paste0(out, sfx, ".index"), warn = FALSE)
    expect_true(any(grepl("^Have completed the analysis of chunk .* bytes [0-9]+$", idx)))
  }
})

test_that("the multi-model scan threads and stays byte-identical", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  man <- mm_manifest(list(binary = "binary_sparseGRM_a", quant = "quant_sparseGRM"))
  out <- file.path(tempfile("step2_mmt_"), "multi.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- mm_run(c(MM_PLINK, MM_COMMON, "--nThreads=4",
                 paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man)), out)
  lg <- mm_log(st)
  # Without this the run could have fallen back to serial and proved nothing.
  expect_true(any(grepl("Marker scan will use 4 in-process threads", lg)))

  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "plink_marker_binary.txt"),
    paste0(out, "_binary")
  )
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "plink_marker_quantitative.txt"),
    paste0(out, "_quant")
  )
})

test_that("a ONE-row manifest is the unsuffixed legacy run", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # The negative control for the whole feature: same machinery, one model,
  # output named and byte-identical to the --GMMATmodelFile run, and NO
  # multi-model log line.
  man <- mm_manifest(list(binary = "binary_sparseGRM_a"))
  out <- file.path(tempfile("step2_mm1_"), "one.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- mm_run(c(MM_PLINK, MM_COMMON,
                 paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man)), out)
  expect_true(file.exists(out))
  expect_false(file.exists(paste0(out, "_binary")))
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "plink_marker_binary.txt"), out
  )
  expect_false(any(grepl("Marker scan will test", mm_log(st))))
})

test_that("two DIFFERENT sample sets get two readers and two decodes", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # No two committed models disagree about sampleID, so one is manufactured:
  # binary_sparseGRM_a with its last 10 samples dropped.
  sub <- local({
    e <- new.env()
    load(file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda"), envir = e)
    m <- e$modglmm
    n <- length(m$sampleID)
    keep <- seq_len(n - 10L)
    # Every field that is indexed by SAMPLE. Taken from the object rather than
    # from a list written here, so a model that grows a field does not silently
    # keep 1000 rows of it beside 990 of everything else -- which is exactly the
    # failure this loop replaced ($Y was missed by hand).
    for (nm in names(m)) {
      x <- m[[nm]]
      if (is.matrix(x)) {
        if (nrow(x) == n) m[[nm]] <- x[keep, , drop = FALSE]
      } else if (is.atomic(x) && length(x) == n) {
        m[[nm]] <- x[keep]
      }
    }
    # obj.noK MUST be rebuilt. ReadModel recomputes it only on the
    # --subSampleFile path (R/readInGLMM.R); the ordinary path trusts what is in
    # the .rda, so a stored 1000-sample obj.noK against a 990-sample model is an
    # "incompatible matrix dimensions: 990x1 and 1000x1" abort in C++. mu and
    # mu2 are not stored either -- ReadModel derives them from fitted.values --
    # so they are derived the same way here.
    mu <- as.vector(m$fitted.values)
    mu2 <- mu * (1 - mu)   # binary; this model is binary_sparseGRM_a
    m$obj.noK <- SAIGE:::ScoreTest_NULL_Model(mu, mu2, m$y, m$X)
    modglmm <- m
    f <- tempfile(fileext = ".rda")
    save(modglmm, file = f)
    f
  })
  vr <- file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt")
  man <- tempfile(fileext = ".txt")
  qp <- mm_model_paths("quant_sparseGRM")
  writeLines(c(paste("binsub", sub, vr), paste("quant", qp[1], qp[2])), man)

  out <- file.path(tempfile("step2_mm2_"), "two.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- mm_run(c(MM_PLINK, MM_COMMON,
                 paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man)), out)
  lg <- mm_log(st)
  expect_true(any(grepl("2 null model\\(s\\) in 2 sample-set group\\(s\\)", lg)))
  expect_true(any(grepl("Marker scan will test 2 models in 2 sample-set group", lg)))

  # The quantitative model is untouched by the other group, so its rows are
  # still golden P32 exactly.
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "plink_marker_quantitative.txt"),
    paste0(out, "_quant")
  )

  # And the two groups really did read different sample sets: the subset model
  # sees 990 samples, the quantitative one 1000. If both had shared a reader
  # this could not hold.
  binsub <- read.table(paste0(out, "_binsub"), header = TRUE, sep = "\t")
  quant <- read.table(paste0(out, "_quant"), header = TRUE, sep = "\t")
  expect_true(all(binsub$N_case + binsub$N_ctrl == 990))
  expect_true(all(quant$N == 1000))

  # The subset model's rows must equal the SAME model run on its own.
  outSolo <- file.path(tempfile("step2_mm2s_"), "solo.txt")
  dir.create(dirname(outSolo), recursive = TRUE)
  mm_run(c(MM_PLINK, MM_COMMON, paste0("--GMMATmodelFile=", sub),
           paste0("--varianceRatioFile=", vr)), outSolo)
  expect_output_matches_golden(outSolo, paste0(out, "_binsub"))
})

test_that("the multi-model scan works through the PGEN reader", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  man <- mm_manifest(list(binary = "binary_sparseGRM_a", quant = "quant_sparseGRM"))
  out <- file.path(tempfile("step2_mmp_"), "multi.txt")
  dir.create(dirname(out), recursive = TRUE)
  mm_run(c(MM_PGEN, MM_COMMON,
           paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man)), out)

  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "pgen_marker_binary.txt"),
    paste0(out, "_binary")
  )
  # No plain quantitative PGEN golden exists, so this one is reconciled against
  # its own single-model run rather than against a committed file.
  outSolo <- file.path(tempfile("step2_mmps_"), "solo.txt")
  dir.create(dirname(outSolo), recursive = TRUE)
  qp <- mm_model_paths("quant_sparseGRM")
  mm_run(c(MM_PGEN, MM_COMMON, paste0("--GMMATmodelFile=", qp[1]),
           paste0("--varianceRatioFile=", qp[2])), outSolo)
  expect_output_matches_golden(outSolo, paste0(out, "_quant"))
})

test_that("conditional analysis composes with a manifest (the LEDGER #40 pin)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # P12's flags. The conditioning factors are assigned per model through that
  # model's own reader, so the binary model's conditional columns must still be
  # exactly golden P12 -- which is the golden that pins the uninitialised-read
  # fix, i.e. the one place a per-model mix-up would show up as garbage.
  man <- mm_manifest(list(binary = "binary_sparseGRM_a", quant = "quant_sparseGRM"))
  out <- file.path(tempfile("step2_mmc_"), "cond.txt")
  dir.create(dirname(out), recursive = TRUE)
  mm_run(c(MM_PLINK, MM_COMMON, "--condition=rs13,rs79",
           paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man)), out)
  expect_output_matches_golden(
    file.path(saige_goldens, "step2", "plink_marker_cond_binary.txt"),
    paste0(out, "_binary")
  )
  # Golden P4 is the quantitative conditional run, but at --minMAC=20, so one
  # command cannot reproduce both it and P12. The quantitative half is
  # reconciled against its own single-model run at these flags instead.
  outSolo <- file.path(tempfile("step2_mmcs_"), "solo.txt")
  dir.create(dirname(outSolo), recursive = TRUE)
  qp <- mm_model_paths("quant_sparseGRM")
  mm_run(c(MM_PLINK, MM_COMMON, "--condition=rs13,rs79",
           paste0("--GMMATmodelFile=", qp[1]),
           paste0("--varianceRatioFile=", qp[2])), outSolo)
  expect_output_matches_golden(outSolo, paste0(out, "_quant"))
})

test_that("the CLI refuses the two combinations that cannot work", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  man <- mm_manifest(list(a = "binary_sparseGRM_a", b = "quant_sparseGRM"))
  p <- mm_model_paths("binary_sparseGRM_a")
  out <- file.path(tempfile("step2_mmx_"), "x.txt")
  dir.create(dirname(out), recursive = TRUE)

  st <- run_saige_cli("step2_SPAtests.R", c(
    MM_PLINK, MM_COMMON,
    paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man),
    paste0("--GMMATmodelFile=", p[1]),
    paste0("--SAIGEOutputFile=", out)
  ))
  expect_false(identical(as.integer(st), 0L))
  expect_true(any(grepl("mutually exclusive", mm_log(st))))

  st2 <- run_saige_cli("step2_SPAtests.R", c(
    MM_PLINK, MM_COMMON,
    paste0("--GMMATmodel_varianceRatio_multiTraits_File=", man),
    "--groupFile=./input/group_new_snpid.txt",
    paste0("--SAIGEOutputFile=", out)
  ))
  expect_false(identical(as.integer(st2), 0L))
  expect_true(any(grepl("single-variant tests only", mm_log(st2))))
})
