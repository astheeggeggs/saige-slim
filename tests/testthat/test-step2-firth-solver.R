# --firth_solver: the approx Firth fit's hat values, closed form (default since
# 2026-09-06) or QR (SAIGE 1.x). Both are the same penalised Newton iteration
# with the same step control and stopping rule; they differ in rounding only
# (<= 8e-13 relative on BETA/SE at n = 50,000, where the two printed byte-identical
# tables on all 70,709 rows). The golden P31 (plink_marker_firth.txt) is produced
# by the default, so it pins `closed`; this file pins that `qr` still agrees with
# it -- as a PARSED comparison at 1e-10, never bytewise: two algebraically equal
# paths are 1 ULP apart and whether that reaches the printed digits is the
# platform's business (CLAUDE.md, "Never assert bit-identity between two
# different code paths"). A warm start at the null intercept was benchmarked the
# same day and rejected -- 118 of 7,891 fits stopped converging -- so there is
# no third value.
#
# Tier: standard (runs the CLI).

FS_ARGS <- function(out, extra = character(0)) {
  c(
    "--bedFile=./input/genotype_100markers.bed",
    "--bimFile=./input/genotype_100markers.bim",
    "--famFile=./input/genotype_100markers.fam",
    "--AlleleOrder=alt-first",
    "--chrom=1", "--minMAF=0", "--minMAC=1", "--LOCO=FALSE",
    "--is_output_moreDetails=TRUE",
    paste0("--GMMATmodelFile=", file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")),
    paste0("--varianceRatioFile=",
           file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt")),
    "--is_Firth_beta=TRUE", "--pCutoffforFirth=0.05",
    extra,
    paste0("--SAIGEOutputFile=", out)
  )
}

fs_run <- function(extra) {
  out <- file.path(tempfile("step2_firth_solver_"), "marker.txt")
  dir.create(dirname(out), recursive = TRUE)
  st <- run_saige_cli("step2_SPAtests.R", FS_ARGS(out, extra))
  list(status = st, out = out)
}

test_that("--firth_solver=qr reproduces the closed-form golden P31 to 1e-10", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- fs_run("--firth_solver=qr")
  expect_cli_ok(r$status)
  golden <- file.path(saige_goldens, "step2", "plink_marker_firth.txt")
  eq <- tables_equivalence(golden, r$out, tol = 1e-10)
  expect_true(eq$ok, info = paste(eq$problems, collapse = "\n"))
  # The comparison must have had Firth rows to compare: 5 of 70 qualify at 0.05.
  g <- read_saige_table(golden)
  expect_equal(sum(suppressWarnings(as.numeric(g$p.value)) <= 0.05, na.rm = TRUE), 5L)
})

test_that("an unknown --firth_solver is refused before any marker is tested", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  r <- fs_run("--firth_solver=newton")
  expect_false(identical(as.integer(r$status), 0L))
  lg <- readLines(attr(r$status, "log"), warn = FALSE)
  expect_true(any(grepl("firth_solver must be 'qr' or 'closed'", lg, fixed = TRUE)))
  expect_false(file.exists(r$out))
})
