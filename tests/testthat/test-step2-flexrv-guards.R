# FlexRV (F3): the CLI guards -- every way a FlexRV run can be mis-specified
# is a named error or a loud warning, never a silent override (LEDGER #30's
# shape). Split from test-step2-flexrv.R so the two files can run in parallel
# workers; scaffolding in helper-step2-flexrv.R. Tier: standard.

test_that("a FlexRV run needs the group file's score line, under the given name", {
  # no score line at all
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                "--groupFile=./input/group_new_snpid.txt", FX_FLEX, FX_BURDEN), "noscore")
  expect_false(identical(as.integer(r$status), 0L))
  expect_true(fx_log_has(r, "has no score line"))
  # a score line under another name
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, "--flexRV_score=CADD", "--flexRV_maxMAF=0.01",
                "--annotation_in_groupTest=lof:missense", FX_BURDEN), "wrongscore")
  expect_false(identical(as.integer(r$status), 0L))
  expect_true(fx_log_has(r, "score line is named 'REVEL'"))
})

test_that("one annotation group and one max MAF per FlexRV run", {
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, "--flexRV_score=REVEL", "--flexRV_maxMAF=0.01",
                "--annotation_in_groupTest=lof,lof:missense", FX_BURDEN), "twoanno")
  expect_false(identical(as.integer(r$status), 0L))
  expect_true(fx_log_has(r, "ONE annotation group"))
  # an explicit --maxMAF_in_groupTest that disagrees is an error naming both flags
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN,
                "--maxMAF_in_groupTest=0.001"), "maxmafclash")
  expect_false(identical(as.integer(r$status), 0L))
  expect_true(fx_log_has(r, "--flexRV_maxMAF (0.01); --maxMAF_in_groupTest=0.001"))
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN,
                "--maxMAF_in_groupTest=0.001,0.01"), "maxmafclash2")
  expect_false(identical(as.integer(r$status), 0L))
  expect_true(fx_log_has(r, "conflicts with it"))
  # --maxMAC_in_groupTest cannot add cells either
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN,
                "--maxMAC_in_groupTest=5"), "maxmac")
  expect_false(identical(as.integer(r$status), 0L))
  expect_true(fx_log_has(r, "--maxMAC_in_groupTest cannot add cells"))
})

test_that("an explicit --maxMAF_in_groupTest equal to --flexRV_maxMAF changes nothing", {
  a <- fx_burden_bin()
  fx_ok(a)
  b <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN,
                "--maxMAF_in_groupTest=0.01"), "maxmafsame")
  fx_ok(b)
  expect_identical(readLines(b$path), readLines(a$path))
})

test_that("a LoF label that matches nothing in the file is an error naming the labels present", {
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN,
                "--flexRV_lofAnno=pLoF,lof"), "lofanno")
  expect_false(identical(as.integer(r$status), 0L))
  expect_true(fx_log_has(r, "'pLoF'"))
  expect_true(fx_log_has(r, "match no annotation in the group file"))
  expect_true(fx_log_has(r, "'missense'"))
  # dropping the lof transform removes the need for the label entirely
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN,
                "--flexRV_lofAnno=pLoF", "--flexRV_scoreTransforms=unw;sbeta:2",
                "--flexRV_MAFTransforms=unw"), "lofanno-unused")
  fx_ok(r)
  expect_identical(sum(r$region$Group != "Cauchy"), 4L)   # 2 sets x 2 genes
})

test_that("--condition needs --weights_for_condition on a FlexRV run", {
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN,
                "--condition=rs30,rs79"), "cond-noweights")
  expect_false(identical(as.integer(r$status), 0L))
  expect_true(fx_log_has(r, "--flexRV_score with --condition needs --weights_for_condition"))
})

test_that("a malformed transform spec is rejected with the flag named", {
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN,
                "--flexRV_scoreTransforms=unw;cubic:0.5"), "badspec")
  expect_false(identical(as.integer(r$status), 0L))
  expect_true(fx_log_has(r, "--flexRV_scoreTransforms: unknown transform 'cubic'"))
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, FX_BURDEN,
                "--flexRV_MAFTransforms=unw;mafbeta:1"), "badspec2")
  expect_false(identical(as.integer(r$status), 0L))
  expect_true(fx_log_has(r, "--flexRV_MAFTransforms: 'mafbeta' parameters must be"))
})

test_that("--r.corr=0 warns but runs; a non-paper max MAF warns and renames the beta sets", {
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, FX_FLEX, "--is_fastTest=FALSE",
                "--flexRV_scoreTransforms=unw;sbeta:2", "--flexRV_MAFTransforms=unw;mac:5"),
              "skato-warn")
  fx_ok(r)
  expect_true(fx_log_has(r, "FlexRV is a burden method; --r.corr=0 runs SKAT-O"))
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, "--flexRV_score=REVEL", "--flexRV_maxMAF=0.02",
                "--annotation_in_groupTest=lof:missense", FX_BURDEN,
                "--flexRV_scoreTransforms=unw", "--flexRV_MAFTransforms=unw;mafbeta:2;mac:5"),
              "maxmaf02")
  fx_ok(r)
  expect_true(fx_log_has(r, "--flexRV_maxMAF=0.02 is not the paper's 0.001"))
  sets <- unique(r$region$Weight[r$region$Group != "Cauchy"])
  expect_identical(sets, c("REVEL_unw_unw", "REVEL_unw_mafbeta2@0.02", "REVEL_unw_MAC5"))
})

test_that("a constant score, and a 0/1 score, are warned about per region", {
  fixture <- readLines(file.path(saige_extdata, "input", "group_new_snpid_flexrv.txt"))
  isScore <- grepl(" score:REVEL ", fixture, fixed = TRUE)
  const <- fixture
  const[isScore] <- vapply(fixture[isScore], function(l) {
    p <- strsplit(l, " ")[[1]]
    paste(c(p[1:2], rep("1", length(p) - 2)), collapse = " ")
  }, "")
  cf <- file.path(saige_extdata, "input", "tmp_flexrv_const_score.txt")
  writeLines(const, cf)
  on.exit(unlink(cf), add = TRUE)
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                "--groupFile=./input/tmp_flexrv_const_score.txt", FX_FLEX, FX_BURDEN,
                "--flexRV_scoreTransforms=unw;sbeta:2", "--flexRV_MAFTransforms=unw"), "const")
  fx_ok(r)
  expect_true(fx_log_has(r, "Region GENE1: the score 'REVEL' is constant (1)"))
  expect_true(fx_log_has(r, "Region GENE2: the score 'REVEL' is constant (1)"))
  # and the degeneracy is real: unw and sbeta2 rows carry the same statistics
  reg <- r$region[r$region$Group != "Cauchy", ]
  expect_identical(reg$Pvalue_Burden[reg$Weight == "REVEL_unw_unw"],
                   reg$Pvalue_Burden[reg$Weight == "REVEL_sbeta2_unw"])
  bin01 <- fixture
  bin01[isScore] <- vapply(fixture[isScore], function(l) {
    p <- strsplit(l, " ")[[1]]
    v <- as.numeric(p[-(1:2)])
    paste(c(p[1:2], ifelse(v >= 0.5, "1", "0")), collapse = " ")
  }, "")
  bf <- file.path(saige_extdata, "input", "tmp_flexrv_binary_score.txt")
  writeLines(bin01, bf)
  on.exit(unlink(bf), add = TRUE)
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                "--groupFile=./input/tmp_flexrv_binary_score.txt", FX_FLEX, FX_BURDEN,
                "--flexRV_scoreTransforms=unw;sbeta:2;cbrt:0.5", "--flexRV_MAFTransforms=unw"),
              "bin01")
  fx_ok(r)
  expect_true(fx_log_has(r, "takes only the values 0 and 1"))
  reg <- r$region[r$region$Group != "Cauchy", ]
  expect_identical(reg$Pvalue_Burden[reg$Weight == "REVEL_cbrt0.5_unw"],
                   reg$Pvalue_Burden[reg$Weight == "REVEL_sbeta2_unw"])
})
