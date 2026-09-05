# LEDGER #57, the open half, end to end: a SKAT kernel failure on a cell keeps
# its NA rows (result files byte-identical to before), is counted, is named
# with its message in <SAIGEOutputFile>.skatFailures.txt, and is summarised on
# stdout at the end of the run.
#
# No bundled input makes the kernel fail (an all-zero custom weight line was
# tried: it changes the answer but SKAT still returns p-values), so the failure
# is PLANTED: a child process replaces SKAT:::Met_SKAT_Get_Pvalue with a
# wrapper that errors on every cell with an ODD marker count and calls the
# original otherwise. get_SKAT_pvalue looks the kernel up by `:::` at call
# time, so the real try/record path is what runs. The plant is a function of
# the cell's INPUT, not of call order, which is what lets the serial and
# 4-worker (fork pool, R/SkatoParallel.R) arms be compared. Child process for
# the usual reason: the C++ singletons are one dataset per process.
#
# Tier: standard (four step-2 runs).

skatfail_run <- function(name, patch, nThreads, condition = "") {
  script <- tempfile(fileext = ".R")
  writeLines(c(
    'args <- commandArgs(trailingOnly = TRUE)',
    'suppressPackageStartupMessages(library(SAIGE))',
    'if (args[2] == "1") {',
    '  orig <- get("Met_SKAT_Get_Pvalue", envir = asNamespace("SKAT"))',
    '  planted <- function(Score, Phi, r.corr, method, Score.Resampling = NULL, ...) {',
    '    if (length(Score) %% 2 == 1) stop("planted kernel failure on a ", length(Score), "-marker cell")',
    '    orig(Score = Score, Phi = Phi, r.corr = r.corr, method = method,',
    '         Score.Resampling = Score.Resampling, ...)',
    '  }',
    '  assignInNamespace("Met_SKAT_Get_Pvalue", planted, ns = "SKAT")',
    '}',
    'extra <- if (nzchar(args[4])) list(condition = args[4]) else list()',
    'do.call(SPAGMMATtest, c(list(',
    '  bedFile = "./input/genotype_100markers.bed",',
    '  bimFile = "./input/genotype_100markers.bim",',
    '  famFile = "./input/genotype_100markers.fam",',
    '  AlleleOrder = "alt-first",',
    '  GMMATmodelFile = args[5], varianceRatioFile = args[6],',
    '  sparseGRMFile = args[7], sparseGRMSampleIDFile = args[8],',
    '  SAIGEOutputFile = args[1],',
    '  chrom = "1", min_MAF = 0, min_MAC = 0.5, LOCO = FALSE,',
    '  groupFile = "./input/group_new_snpid.txt",',
    '  annotation_in_groupTest = c("lof", "missense;lof",',
    '                              "missense;lof;synonymous"),',
    '  maxMAF_in_groupTest = c(0.0001, 0.001, 0.01),',
    '  is_fastTest = TRUE, nThreads = as.integer(args[3])), extra))'
  ), script)
  s1 <- file.path(saige_extdata, "output", "example_binary_sparseGRM")
  sg <- file.path(saige_extdata, "output",
    "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
  out <- tempfile(pattern = name, fileext = ".txt")
  log <- tempfile(fileext = ".log")
  owd <- setwd(saige_extdata); on.exit(setwd(owd), add = TRUE)
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c(script, out, patch, nThreads, shQuote(condition),
      paste0(s1, ".rda"), paste0(s1, ".varianceRatio.txt"),
      sg, paste0(sg, ".sampleIDs.txt")),
    stdout = log, stderr = log))
  expect_cli_ok(st)
  lines <- readLines(log, warn = FALSE)
  res <- utils::read.table(out, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                           na.strings = "NA", check.names = FALSE)
  side <- paste0(out, ".skatFailures.txt")
  list(path = out, log = lines, res = res, side = side,
       summary = grep("^SKAT failures: ", lines, value = TRUE),
       by_msg = grep("^  [0-9]+ x ", lines, value = TRUE),
       sidecar = if (file.exists(side)) {
         utils::read.table(side, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                           na.strings = "NA")
       } else NULL)
}

na_cells <- function(res) sum(res$Group != "Cauchy" & is.na(res$Pvalue))

test_that("NEGATIVE CONTROL: an unpatched run reports 0, writes no sidecar, and is golden P2", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- skatfail_run("neg", 0, 1)
  expect_identical(r$summary, "SKAT failures: 0 cells returned NA from the SKAT kernel.")
  expect_null(r$sidecar)
  expect_identical(na_cells(r$res), 0L)
  expect_output_matches_golden(file.path(saige_goldens, "step2", "plink_group_fastTest.txt"), r$path)
})

test_that("a planted kernel failure keeps its NA rows, and every one is counted and named", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- skatfail_run("pos", 1, 1)
  n_na <- na_cells(r$res)
  expect_gt(n_na, 0)                               # the plant actually fired
  expect_false(is.null(r$sidecar))
  expect_identical(nrow(r$sidecar), n_na)          # one sidecar row per NA cell
  expect_match(r$summary, paste0("^SKAT failures: ", n_na, " cell\\(s\\) in 2 region\\(s\\)"))
  expect_true(all(r$sidecar$Test == "marginal"))
  expect_true(all(grepl("^planted kernel failure on a [0-9]+-marker cell$", r$sidecar$Error)))
  # the sidecar names exactly the NA cells, by (Region, Group, max_MAF)
  na_rows <- r$res[r$res$Group != "Cauchy" & is.na(r$res$Pvalue), ]
  key <- function(d) sort(paste(d$Region, d$Group, as.character(d$max_MAF)))
  expect_identical(key(r$sidecar), key(na_rows))
  # the by-message tally sums to the count and is what the sidecar says
  counts <- as.integer(sub("^  ([0-9]+) x .*$", "\\1", r$by_msg))
  expect_identical(sum(counts), n_na)
  expect_identical(sort(counts, decreasing = TRUE), counts)
  # unaffected cells and the region's Cauchy row are untouched: finite p-values
  expect_true(all(is.finite(r$res$Pvalue[!is.na(r$res$Pvalue)])))
  expect_false(any(is.na(r$res$Pvalue[r$res$Group == "Cauchy"])))
})

test_that("the failure field survives the SKAT-O fork pool: 4 workers match serial", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  skip_on_os("windows")
  serial <- skatfail_run("pos1", 1, 1)
  par <- skatfail_run("pos4", 1, 4)
  expect_true(any(grepl("use 4 workers", par$log, fixed = TRUE)))   # it really pooled
  expect_identical(readLines(par$path), readLines(serial$path))     # same NA rows, same bytes
  # the summary line embeds the sidecar path, which differs per run
  strip_path <- function(s) sub(" in '[^']*'", " in '<path>'", s)
  expect_identical(strip_path(par$summary), strip_path(serial$summary))
  expect_identical(par$by_msg, serial$by_msg)
  o <- function(d) d[order(d$Region, d$Group, d$max_MAF, d$Test), ]
  expect_identical(o(par$sidecar), o(serial$sidecar))
})

test_that("the conditional call site is counted too, labelled 'conditional'", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- skatfail_run("cond", 1, 1, condition = "rs30,rs79")
  expect_false(is.null(r$sidecar))
  tt <- table(r$sidecar$Test)
  expect_identical(sort(names(tt)), c("conditional", "marginal"))
  # the same cells fail both ways (the plant keys on the cell's marker count)
  expect_identical(tt[["conditional"]], tt[["marginal"]])
  expect_identical(tt[["marginal"]], na_cells(r$res))
  cond_na <- sum(r$res$Group != "Cauchy" & is.na(r$res$Pvalue_cond))
  expect_identical(cond_na, tt[["conditional"]])
})
