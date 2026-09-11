# FlexRV: knobs and paths that must not move the answer, each asserted to
# have ARRIVED (a knob that fails to reach the run passes as green otherwise,
# the LEDGER #30 shape). Byte identity is claimed only for ONE path re-run
# with an inert knob; a relabelled fixture is compared after the label is
# mapped back. Scaffolding in helper-step2-flexrv.R. Tier: standard.

test_that("--flexRV_lofAnno: a fixture whose LoF label is 'pLoF' reproduces the default run", {
  fx_skip("standard")
  base <- fx_burden_bin()
  fx_ok(base)
  lines <- fx_fixture_lines()
  isAnno <- grepl(" anno ", lines, fixed = TRUE)
  lines[isAnno] <- gsub(" lof", " pLoF", lines[isAnno], fixed = TRUE)
  expect_false(any(grepl(" lof", lines[isAnno], fixed = TRUE)))
  gf <- tempfile("flexrv_plof_", fileext = ".txt")
  writeLines(lines, gf)
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                paste0("--groupFile=", gf), "--flexRV_score=REVEL", "--flexRV_maxMAF=0.01",
                "--annotation_in_groupTest=pLoF:missense", "--flexRV_lofAnno=pLoF", FX_BURDEN),
              "plof")
  fx_ok(r)
  d <- r$region
  expect_true(all(d$Group[d$Group != "Cauchy"] == "pLoF;missense"))
  d$Group[d$Group == "pLoF;missense"] <- "lof;missense"
  expect_identical(d, base$region)
})

test_that("--flexRV_lofAnno=lof,missense makes the lof mask all ones: lof_* rows are unw_* rows", {
  fx_skip("standard")
  base <- fx_burden_bin()
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN,
                "--flexRV_lofAnno=lof,missense"), "lofboth")
  fx_ok(r)
  rows <- fx_rows(r)
  strip <- function(d) {
    d$Weight <- NULL
    rownames(d) <- NULL
    d
  }
  mafs <- unique(sub("^REVEL_unw_", "", grep("^REVEL_unw_", rows$Weight, value = TRUE)))
  expect_identical(length(mafs), 12L)
  for (mf in mafs) {
    expect_identical(strip(rows[rows$Weight == paste0("REVEL_lof_", mf), ]),
                     strip(rows[rows$Weight == paste0("REVEL_unw_", mf), ]), label = mf)
  }
  # the other 14 score transforms do not read the LoF label: their rows are the
  # default run's; and the default run's lof_* rows are NOT the unw_* rows
  b <- fx_rows(base)
  keep <- !grepl("^REVEL_lof_", rows$Weight)
  expect_identical(rows[keep, ], b[keep, ])
  expect_false(identical(strip(b[b$Weight == "REVEL_lof_unw", ]),
                         strip(b[b$Weight == "REVEL_unw_unw", ])))
})

test_that("--markers_per_chunk_in_groupTest=10 arrives and is byte-inert on the 192-set BURDEN", {
  # GENE1 has 50 markers, so 10 per chunk is five chunks plus the ultra-rare
  # one; the per-marker transform weight and the per-set collapse are both
  # chunk-local, so the bytes cannot move (finding E0)
  fx_skip("standard")
  base <- fx_burden_bin()
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN,
                "--markers_per_chunk_in_groupTest=10"), "chunk10")
  fx_ok(r)
  at <- which(r$log == "$markers_per_chunk_in_groupTest")
  expect_length(at, 1L)
  expect_identical(r$log[at + 1], "[1] 10")
  expect_identical(readLines(r$path), readLines(base$path))
})

test_that("SKAT-O transform sets: --nThreads=4 is byte-identical to one thread, and says so", {
  # the transform descriptor through the forked SKAT-O pool and the
  # identical-cell memo (R/SkatoParallel.R); the served count is a
  # MEASUREMENT (printed-equal is not bit-equal), so only the arithmetic of
  # the line is asserted
  fx_skip("standard")
  one <- fx_skato()
  fx_ok(one)
  four <- fx_skato("--nThreads=4", "skato-threads4")
  fx_ok(four)
  expect_true(fx_log_has(four, "Region marker loop will use 4 in-process threads"))
  expect_true(fx_log_has(four, "SKAT-O post-processing both use 4 threads/workers"))
  for (res in list(one, four)) {
    l <- grep("^SKAT-O cells: ", res$log, value = TRUE)
    expect_length(l, 1L)
    n <- as.integer(regmatches(l, gregexpr("[0-9]+", l))[[1]])   # asked, evaluated, served
    expect_identical(n[1], 8L)                                    # 2 genes x 4 sets
    expect_identical(n[2] + n[3], n[1])
  }
  expect_identical(readLines(four$path), readLines(one$path))
  expect_identical(readLines(paste0(four$path, ".singleAssoc.txt")),
                   readLines(paste0(one$path, ".singleAssoc.txt")))
  expect_false(any(fx_rows(one)$Pvalue == "NA"))
})
