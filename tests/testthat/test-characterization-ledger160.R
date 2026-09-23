# LEDGER #160 (FIXED): a region with exactly ONE marker no longer makes R and
# C++ disagree on how many weight sets it composed.
#
# R hands C++ the region's weights through weight_matrix_for_CPP. A group file
# with no `weight` line sets the sentinel `WEIGHT <- c(0)`, which that function
# shaped as a 1 x 1 matrix; a genuine one-set weight line on a one-marker
# region is a q0 x 1 = 1 x 1 matrix too. C++ decides `hasCustomWeights` on
# `t_weight.n_rows == q0` ALONE -- #146(b) removed the `!t_weight.is_zero()`
# half that had been telling the two apart -- so at q0 == 1 it read the
# sentinel as ONE custom set. R had composed only the transform sets, C++
# composed one more, and #146(a)'s assert stopped the run. Every region with
# two or more markers was unaffected, which is why no bundled fixture and no
# golden reached it: it needs a gene whose every-but-one variant is filtered
# out, which depends on the cohort.
#
# The sentinel is now a 0 x 0 matrix: no columns, so `q_weight_custom` is 0 at
# every q0. A genuine q0 x K weight matrix is untouched, including at q0 == 1.
#
# Fixture derived at run time: the bundled FlexRV group file plus GENE3, one
# variant (rs4, a real marker at the default collapse cutoff). Tier: standard.

l160_fixture <- function(dir) {
  src <- file.path(get("saige_extdata"), "input", "group_new_snpid_flexrv.txt")
  path <- file.path(dir, "group_one_marker.txt")
  writeLines(c(readLines(src),
               "GENE3 var rs4",
               "GENE3 anno missense",
               "GENE3 score:REVEL 0.5"), path)
  path
}

test_that("a one-marker region composes its sets and does not stop (LEDGER #160)", {
  fx_skip("standard")
  dir <- tempfile("step2_l160_")
  dir.create(dir, recursive = TRUE)
  gf <- l160_fixture(dir)

  common <- FX_COMMON[!grepl("--groupFile=", FX_COMMON, fixed = TRUE)]
  res <- fx_run(c(FX_G, FX_BIN, FX_SG, common, paste0("--groupFile=", gf),
                  FX_FLEX, FX_BURDEN), "l160")

  # Before the fix this was a stop naming "C++ composed <K+1> (custom 1, ...)".
  fx_ok(res)
  expect_false(any(grepl("mainRegionInCPP: R composed", res$log)))

  reg <- res$region
  expect_true("GENE3" %in% reg$Region)

  # The one-marker region's sets are TRANSFORM sets. A phantom custom set would
  # show up as a "W<r>" label, which is the fallback C++ uses when it has more
  # sets than R gave it names for.
  g3 <- reg[reg$Region == "GENE3" & reg$Group != "Cauchy", ]
  expect_gt(nrow(g3), 0)
  expect_true(all(grepl("^REVEL_", g3$Weight)))
  expect_false(any(grepl("^W[0-9]+$", g3$Weight)))

  # The fix must be inert on every region that already worked: GENE1 and GENE2
  # come from the untouched fixture and must match that run cell for cell.
  base <- fx_burden_bin()
  fx_ok(base)
  key <- function(d) {
    d <- d[d$Region %in% c("GENE1", "GENE2"), ]
    d[order(d$Region, d$Group, d$max_MAF, d$Weight), ]
  }
  expect_equal(key(reg), key(base$region), ignore_attr = TRUE)
})
