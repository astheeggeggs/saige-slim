# FlexRV: a group file may carry SEVERAL consecutive score lines (2026-10-27);
# --flexRV_score names the one a run uses and the others change no byte. Why:
# the October BRaVa release carried three FlexRV files per chromosome that
# differed only in their score line (two predictors and the random-score null),
# the identical var/anno lines duplicated three times by the one-line rule.
# Fixture group_new_snpid_flexrv_multi.txt is group_new_snpid_flexrv.txt with a
# second line, score:SHUF, the REVEL values reversed within each region.
#   (1) REVEL selected from the two-score file reproduces golden P58, with the
#       score line first or second;
#   (2) SHUF selected gives a different answer under SHUF-named sets;
#   (3) a name the file does not carry is refused, naming both lines;
#   (4) with no FlexRV flag the extra line is inert: a SKAT-O run on the
#       two-score file equals the same run on the one-score file;
#   (5) the parser's layout rules: distinct names, together after anno.
# Scaffolding in helper-step2-flexrv.R. Tier: standard.

fx_multi_file <- file.path(saige_extdata, "input", "group_new_snpid_flexrv_multi.txt")
FX_MULTI <- sub("group_new_snpid_flexrv.txt", "group_new_snpid_flexrv_multi.txt", FX_COMMON,
                fixed = TRUE)
fx_p58 <- file.path(saige_goldens, "step2", "plink_group_flexrv_burden.txt")

# the fixture with its two score lines swapped in every region (SHUF first)
fx_swapped_file <- function() {
  lines <- readLines(fx_multi_file)
  is_score <- grepl("^\\S+ score:", lines)
  idx <- which(is_score)
  stopifnot(length(idx) == 4, all(diff(idx)[c(1, 3)] == 1))
  swapped <- lines
  swapped[idx[c(1, 2, 3, 4)]] <- lines[idx[c(2, 1, 4, 3)]]
  f <- tempfile(fileext = ".txt")
  writeLines(swapped, f)
  f
}

test_that("REVEL selected from the two-score file reproduces golden P58, first or second", {
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_MULTI, FX_FLEX, FX_BURDEN), "multi_revel")
  fx_ok(r)
  expect_output_matches_golden(fx_p58, r$path)
  expect_true(fx_log_has(r, "from score 'REVEL'"))
  args_sw <- c(FX_G, FX_BIN, FX_SG, FX_COMMON[!grepl("--groupFile=", FX_COMMON)],
               paste0("--groupFile=", fx_swapped_file()), FX_FLEX, FX_BURDEN)
  r2 <- fx_run(args_sw, "multi_revel_second")
  fx_ok(r2)
  expect_output_matches_golden(fx_p58, r2$path)
})

test_that("SHUF selected gives SHUF-named sets and a different answer", {
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_MULTI, sub("REVEL", "SHUF", FX_FLEX, fixed = TRUE),
                FX_BURDEN), "multi_shuf")
  fx_ok(r)
  expect_true(fx_log_has(r, "from score 'SHUF'"))
  rows <- fx_rows(r)
  expect_true(all(grepl("^SHUF_", rows$Weight)))
  expect_false(any(grepl("^REVEL_", rows$Weight)))
  g <- fx_read(fx_p58)
  gr <- g[g$Group != "Cauchy", ]
  expect_equal(nrow(rows), nrow(gr))
  # the same cell grid under another score: the unweighted sets agree, the score-driven ones
  # (the sbeta sets, where REVEL's 1s and SHUF's reversed values differ most) do not
  key <- function(d) paste(d$Region, d$Group, d$max_MAF, sub("^[A-Z]+_", "", d$Weight))
  m <- match(key(rows), key(gr))
  expect_false(anyNA(m))
  p_new <- as.numeric(rows$Pvalue_Burden)
  p_old <- as.numeric(gr$Pvalue_Burden[m])
  unw <- grepl("_unw_", rows$Weight)
  expect_equal(p_new[unw], p_old[unw], tolerance = 1e-10)
  expect_true(any(abs(p_new[!unw] - p_old[!unw]) > 1e-6 * pmax(p_old[!unw], 1e-300)))
})

test_that("a score the file does not carry is refused, naming both lines", {
  r <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_MULTI, sub("REVEL", "CADD", FX_FLEX, fixed = TRUE),
                FX_BURDEN), "multi_wrong")
  expect_false(identical(as.integer(r$status), 0L))
  expect_true(fx_log_has(r, "score lines are named 'REVEL', 'SHUF'"))
})

test_that("with no FlexRV flag the extra score line is inert", {
  skato <- c("--annotation_in_groupTest=lof,missense:lof", "--maxMAF_in_groupTest=0.01",
             "--is_fastTest=TRUE")
  a <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_COMMON, skato), "one_score_skato")
  b <- fx_run(c(FX_G, FX_BIN, FX_SG, FX_MULTI, skato), "two_score_skato")
  fx_ok(a)
  fx_ok(b)
  # one path, re-run on an input that differs only in a line nothing reads
  expect_output_matches_golden(a$path, b$path)
})

test_that("the parser takes consecutive score lines with distinct names, after anno only", {
  d <- checkGroupFile(fx_multi_file, detail = TRUE)
  expect_identical(d$scoreName, c("REVEL", "SHUF"))
  lines <- readLines(fx_multi_file)
  dup <- sub("score:SHUF", "score:REVEL", lines, fixed = TRUE)
  f <- tempfile()
  writeLines(dup, f)
  expect_error(checkGroupFile(f), "Duplicate score name 'REVEL'")
  # the second region carrying the lines in the other order is refused
  idx <- which(grepl("^GENE2 score:", lines))
  swapped2 <- lines
  swapped2[idx] <- lines[rev(idx)]
  f2 <- tempfile()
  writeLines(swapped2, f2)
  expect_error(checkGroupFile(f2), "same score lines in the same order")
})
