# LEDGER #172 (2026-09-29, FIXED): a weighted set whose weight ONE ATOM carries
# takes the exact test on its material carriers, whatever its allele count.
# Two mechanisms made the saddlepoint converge to 1e-53..1e-93 on such cells
# where ~1e-3 was right (the variance-ratio rescaling crossing the atom; a
# pinned sample's phantom atom); the exact test enumerates the coin flips and
# ignores the variance ratio, as it already did for tiny unweighted sets.
# Routing: Kish's effective carrier count on the vector the test sees at or
# below --max_MAC_for_ER, with at most 13 material carriers (the enumeration
# regime), for custom and transform sets on both region paths.
#   1. P69 (FlexRV atom fixture): the golden; the run names every routed cell
#      in `.exactByWeight.txt` and prints the line; every routed set is
#      sbeta10/sbeta20 and every routed cell has the atom's effective count.
#   2. SAIGE_EXACT_BY_WEIGHT=0 restores the saddlepoint: the rows that differ
#      from the golden are exactly the routed cells plus the Cauchy rows that
#      combine them, and the saddlepoint's p on the singleton atom is an order
#      of magnitude below the exact test's (the artefact, in miniature).
#   3. One statistic two ways: each routed cell's p equals an R enumeration of
#      the exact test (mid-p: half the observed configuration plus the more
#      extreme ones, SAIGE's convention) over the material carriers read from
#      the saddlepoint dump of the switched-off run; and the singleton atom's
#      routed p equals P1's own exact-test row for that marker.
#   4. P70/P71 (custom-weights atom fixture): both region paths route (sidecar
#      kinds "burden" and "ultra_rare"), goldens match.
skip_if_no_goldens()

common_args <- c("--bedFile=./input/genotype_100markers.bed",
                 "--bimFile=./input/genotype_100markers.bim",
                 "--famFile=./input/genotype_100markers.fam",
                 "--AlleleOrder=alt-first",
                 "--chrom=1", "--minMAF=0", "--minMAC=0.5", "--LOCO=FALSE",
                 "--GMMATmodelFile=./output/example_binary_sparseGRM.rda",
                 "--varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt",
                 paste0("--sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
                        "randomMarkersUsed.sparseGRM.mtx"),
                 paste0("--sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_",
                        "randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"),
                 "--is_fastTest=FALSE")
flexrv_args <- function(out) {
  c(common_args, paste0("--SAIGEOutputFile=", out),
    "--groupFile=./input/group_new_snpid_flexrv_atom.txt",
    "--flexRV_score=REVEL", "--flexRV_maxMAF=0.01",
    "--annotation_in_groupTest=lof:missense", "--r.corr=1")
}
custom_args <- function(out, burden) {
  c(common_args, paste0("--SAIGEOutputFile=", out),
    "--groupFile=./input/group_new_snpid_withWeights_atom.txt",
    "--annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous",
    "--maxMAF_in_groupTest=0.0001,0.001,0.01",
    if (burden) "--r.corr=1" else character())
}
gold <- function(f) file.path(saige_goldens, "step2", f)
routed_line <- function(log) grep("^Exact test by weight: ", log, value = TRUE)
routed_count <- function(log) {
  as.integer(sub("^Exact test by weight: ([0-9]+) .*$", "\\1", routed_line(log)))
}
read_sidecar <- function(path) {
  utils::read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

test_that("P69: the FlexRV atom fixture routes every sbeta10/sbeta20 cell and matches the golden", {
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", flexrv_args(out))
  expect_cli_ok(st)
  log <- readLines(attr(st, "log"), warn = FALSE)
  expect_output_matches_golden(gold("plink_group_flexrv_atom.txt"), out)
  side <- read_sidecar(paste0(out, ".exactByWeight.txt"))
  expect_identical(routed_count(log), nrow(side))
  expect_gt(nrow(side), 0)
  expect_true(all(side$Kind == "burden"))
  expect_true(all(grepl("_sbeta(10|20)_", side$Weight)))
  # the atom: rs44 is a singleton (GENE1), rs56 a doubleton (GENE2); the
  # effective count is the atom's allele count to within the negligible carriers
  expect_true(all(abs(side$n_eff[side$Region == "GENE1"] - 1) < 1e-3))
  expect_true(all(abs(side$n_eff[side$Region == "GENE2"] - 2) < 1e-3))
  expect_true(all(side$MAC > 4))     # above the allele-count cutoff: the MAC rule never routes them
  expect_true(all(side$n_material <= 13))
  gs <- read_sidecar(gold("plink_group_flexrv_atom.txt.exactByWeight.txt"))
  keep <- c("Kind", "Region", "Group", "Weight")
  expect_identical(side[, keep], gs[, keep])
})

test_that("the kill switch restores the saddlepoint on the routed cells only (the artefact)", {
  out <- tempfile(fileext = ".txt")
  st <- run_saige_cli("step2_SPAtests.R", flexrv_args(out),
                      env = c(SAIGE_EXACT_BY_WEIGHT = "0"))
  expect_cli_ok(st)
  log <- readLines(attr(st, "log"), warn = FALSE)
  expect_identical(routed_count(log), 0L)
  expect_false(file.exists(paste0(out, ".exactByWeight.txt")))
  new <- read_saige_table(gold("plink_group_flexrv_atom.txt"))
  old <- read_saige_table(out)
  side <- read_sidecar(gold("plink_group_flexrv_atom.txt.exactByWeight.txt"))
  expect_identical(dim(old), dim(new))
  changed <- which(old$Pvalue_Burden != new$Pvalue_Burden)
  cells <- new$Group != "Cauchy"
  routed_key <- paste(side$Region, side$Weight)
  expect_setequal(paste(new$Region, new$Weight)[changed[cells[changed]]], routed_key)
  expect_true(all(new$Group[changed[!cells[changed]]] == "Cauchy"))
  # every non-routed cell is byte-identical between the two runs
  unrouted_old <- cells & !(paste(old$Region, old$Weight) %in% routed_key)
  unrouted_new <- cells & !(paste(new$Region, new$Weight) %in% routed_key)
  expect_identical(old[unrouted_old, ], new[unrouted_new, ])
  # the artefact in miniature: on the singleton atom the saddlepoint sits an
  # order of magnitude below the exact test (1.8e-3 against 3.2e-2 on this fixture)
  i <- which(new$Region == "GENE1" & new$Weight == "REVEL_sbeta20_unw")
  expect_lt(as.numeric(old$Pvalue_Burden[i]) * 10, as.numeric(new$Pvalue_Burden[i]))
})

# the exact test's mid-p over the material carriers, enumerated
enum_midp <- function(Z, mu, y_obs, eps = 1e-6) {
  k <- length(Z)
  Q_obs <- sum(Z * (y_obs - mu))^2
  p <- 0
  for (cfg in 0:(2^k - 1)) {
    y <- as.integer(intToBits(cfg))[seq_len(k)]
    Q <- sum(Z * (y - mu))^2
    pr <- prod(ifelse(y == 1, mu, 1 - mu))
    if (Q > Q_obs + eps) p <- p + pr else if (Q >= Q_obs - eps) p <- p + 0.5 * pr
  }
  p
}
parse_kv <- function(line) {
  f <- strsplit(sub("^# ", "", line), "\t", fixed = TRUE)[[1]]
  stats::setNames(f[seq(2, length(f), 2)], f[seq(1, length(f), 2)])
}

test_that("one statistic two ways: each routed cell's exact p equals the enumeration", {
  out <- tempfile(fileext = ".txt")
  dump <- tempfile("spadump_")
  dir.create(dump)
  st <- run_saige_cli("step2_SPAtests.R", flexrv_args(out),
                      env = c(SAIGE_EXACT_BY_WEIGHT = "0", SAIGE_SPA_DUMP_DIR = dump,
                              SAIGE_SPA_DUMP_BELOW = "1"))
  expect_cli_ok(st)
  new <- read_saige_table(gold("plink_group_flexrv_atom.txt"))
  side <- read_sidecar(gold("plink_group_flexrv_atom.txt.exactByWeight.txt"))
  old <- read_saige_table(out)
  files <- list.files(dump, pattern = "^spaBelow_", full.names = TRUE)
  expect_gt(length(files), 0)
  beta_d <- vapply(files, function(f) {
    h <- parse_kv(readLines(f, n = 1))
    as.numeric(h[["Tstat"]]) / as.numeric(h[["var1"]])
  }, numeric(1))
  n_checked <- 0L
  for (i in seq_len(nrow(side))) {
    sel_old <- old$Region == side$Region[i] & old$Weight == side$Weight[i] & old$Group != "Cauchy"
    row_old <- old[sel_old, ]
    beta <- as.numeric(row_old$BETA_Burden)
    cand <- which(abs(beta_d - beta) <= 2e-5 * abs(beta))
    expect_gt(length(cand), 0)
    # cells sharing the atom print the same BETA and differ only in negligible
    # carriers: the dump whose material count is the sidecar's is this cell's
    for (j in cand) {
      d <- utils::read.table(files[j], header = TRUE, skip = 2, sep = "\t")
      nz <- d$nz == 1
      zmax <- max(abs(d$g[nz]))
      Z <- if (zmax < 1) d$g / zmax else d$g
      material <- nz & abs(Z) >= 1e-6
      if (sum(material) != side$n_material[i]) next
      p_enum <- enum_midp(Z[material], d$mu[material], d$y[material])
      sel_new <- new$Region == side$Region[i] & new$Weight == side$Weight[i] & new$Group != "Cauchy"
      p_new <- as.numeric(new$Pvalue_Burden[sel_new])
      # the exact test adjusts the configuration probabilities slightly (p1_adj);
      # measured 4e-4..4e-3 relative on these cells, against saddlepoints 18x off
      expect_equal(p_new, p_enum, tolerance = 1e-2)
      n_checked <- n_checked + 1L
      break
    }
  }
  expect_identical(n_checked, nrow(side))
  # the singleton atom IS rs44: the routed sets print P1's exact-test row for it
  # (a different null model file of the same fit family; 6 printed digits)
  p1 <- read_saige_table(gold("plink_marker_binary.txt"))
  p_rs44 <- as.numeric(p1$p.value[p1$MarkerID == "rs44"])
  p_set <- as.numeric(new$Pvalue_Burden[new$Region == "GENE1" & new$Weight == "REVEL_sbeta20_MAC3"])
  expect_equal(p_set, p_rs44, tolerance = 2e-2)
})

test_that("P70/P71: the custom-weights atom fixture routes on both region paths", {
  for (burden in c(TRUE, FALSE)) {
    out <- tempfile(fileext = ".txt")
    st <- run_saige_cli("step2_SPAtests.R", custom_args(out, burden))
    expect_cli_ok(st)
    log <- readLines(attr(st, "log"), warn = FALSE)
    g <- paste0("plink_group_customweight_atom_", if (burden) "burden" else "skato", ".txt")
    expect_output_matches_golden(gold(g), out)
    side <- read_sidecar(paste0(out, ".exactByWeight.txt"))
    expect_identical(routed_count(log), nrow(side))
    expect_gt(nrow(side), 0)
    expect_true(all(side$Kind == if (burden) "burden" else "ultra_rare"))
    expect_true(all(is.na(side$Weight) | side$Weight == ""))   # a plain weight line has no set name
    expect_true(all(side$n_material <= 2))     # the atom alone (rs44) or the doubleton (rs56)
    expect_true(all(side$MAC > 4))
    if (!burden) expect_length(grep("^SKAT failures: 0", log, value = TRUE), 1L)
  }
})
