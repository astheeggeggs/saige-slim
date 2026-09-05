# LEDGER #71 FIXED: the two R-API entry points defaulted the categorical MAC
# bin edges differently (step 1: c(10, 20.5); step 2: c(10.5, 20.5)) while
# both CLI drivers pass "10,20.5". The step-2 default now matches, and the
# roxygen for both says what the code does. Fast tier.

cli_default <- function(script, option) {
  src <- readLines(file.path(saige_extdata, script), warn = FALSE)
  i <- grep(paste0('"--', option, '"'), src, fixed = TRUE)
  expect_length(i, 1)
  blk <- paste(src[i:(i + 6)], collapse = " ")
  sub('.*default = "([^"]*)".*', "\\1", blk)
}

test_that("cateVarRatioMinMACVecExclude defaults agree across step 1, step 2 and both CLIs", {
  s1 <- eval(formals(SAIGE::fitNULLGLMM)$cateVarRatioMinMACVecExclude)
  s2 <- eval(formals(SAIGE::SPAGMMATtest)$cateVarRatioMinMACVecExclude)
  expect_identical(s1, c(10, 20.5))
  expect_identical(s2, s1)
  expect_identical(cli_default("step1_fitNULLGLMM.R", "cateVarRatioMinMACVecExclude"), "10,20.5")
  expect_identical(cli_default("step2_SPAtests.R", "cateVarRatioMinMACVecExclude"), "10,20.5")
})
